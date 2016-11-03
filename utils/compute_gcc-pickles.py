#!/usr/bin/env python

import os
import subprocess
import sys
from joblib import Parallel, delayed

def main():
    input=sys.argv[1]
    output=sys.argv[2]
    consampy="/home/acormier/working_directory/ionplugins/IonNIPT/wisecondor/consam.py"
    gccpy="/home/acormier/working_directory/ionplugins/IonNIPT/wisecondor/gcc.py"
    gccount="/home/acormier/working_directory/ionplugins/IonNIPT/data/hg19.gccount"
    
    bam_lst = []
    for bam in os.listdir(input):
        if bam.endswith(".bam") or bam.endswith(".cram"):
            bam_lst.append(bam)
    Parallel(n_jobs=30)(delayed(getGCCPickles)(bam, input, output, consampy, gccpy, gccount) for bam in bam_lst)
    
def getGCCPickles(bam, input, output, consampy, gccpy, gccount):
    
    name = os.path.basename(bam).split('.')[0]
    bam=os.path.join(input, bam)
    pickle=os.path.join(output, name+'.pickle')
    gcc=os.path.join(output, name+'.gcc')
    
    pickle_cmd = "samtools view {bam} -q 1 | python {consampy} -outfile {pickle}".format(bam = bam, consampy = consampy, pickle = pickle)
    process = subprocess.Popen(pickle_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(stderr)
            
    gcc_cmd = "python {gccpy} {pickle} {gccount} {gcc}".format(gccpy = gccpy, pickle = pickle, gccount = gccount, gcc = gcc)
    process = subprocess.Popen(gcc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(stderr)

if __name__ == '__main__':
    main()