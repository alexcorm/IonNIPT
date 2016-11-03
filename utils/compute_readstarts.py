#!/usr/bin/env python

import os
import subprocess
import sys
from joblib import Parallel, delayed

def main():
    input=sys.argv[1]
    output=sys.argv[2]
    retro="/home/acormier/working_directory/ionplugins/IonNIPT/sanefalcon/retro.py"
    
    for bam in os.listdir(input):
        if bam.endswith(".bam") or bam.endswith(".cram"):
            bam_path=os.path.join(input, bam)
            Parallel(n_jobs=23)(delayed(getReadStarts)(bam_path, id, output, retro) for id in range(1,23))

def getReadStarts(bam, id, output, retro):
    name = os.path.basename(bam).split('.')[0]
    
    cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
    cmd_forward += "'{print $4}' > "
    cmd_forward += "{output}/{name}.{id}.start.fwd".format(output=output,name=name, id=id)
    
    cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
    cmd_reverse += "'{print ($4 + length($10) - 1)}' "
    cmd_reverse += "> {output}/{name}.{id}.start.rev".format(output=output, name=name, id=id)
    
    jobLauncher(cmd_forward)
    jobLauncher(cmd_reverse)
    
def jobLauncher(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise Exception(stderr)

if __name__ == '__main__':
    main()