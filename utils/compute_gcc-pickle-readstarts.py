#!/usr/bin/env python
# -*- coding: utf-8 -*-

#================================================================#

import os
import sys
import re
import json
import subprocess
from joblib import Parallel, delayed

#================================================================#

def main():
    
    #
    NIPT="/home/ionadmin/acormier/dpni/IonNIPT"
    out="/home/ionadmin/acormier/dpni/test"
    
    consampy="%s/wisecondor/consam.py" % NIPT
    gccpy="%s/wisecondor/gcc.py" % NIPT
    gccount="%s/data/hg19.gccount" % NIPT
    retropy="%s/sanefalcon/retro.py" % NIPT
    
    # load run info
    path_run = sys.argv[1] # folder where bam and json files are located
    dpniDict = {}

    json_file = '%s/ion_params_00.json' % path_run
    json_load = json.load(open(json_file))
    runname = json_load['exp_json']['log']['runname'] #.split('-DANNI')[0]
    
    dpniDict[runname] = []
    
    #get sample and barcode name
    for sample, barcode in json_load['experimentAnalysisSettings']['barcodedSamples'].items():
        
        sample_name = sample.replace(' ', '_') # just in case...
        barcode_name = barcode['barcodeSampleInfo'].keys()[0]
        
        dpniDict[runname].append([barcode_name, sample_name])
    
    for run, design in dpniDict.items():
        for sample in design:
            
            barcode_id = sample[0]
            sample_name = sample[1]
            
            name = run + "_" + sample_name
            
            bam = os.path.join(path_run, barcode_id)+"_rawlib.bam"
            readstarts = os.path.join(out, name)
            pickle = os.path.join(out, name)+".pickle"
            gcc = os.path.join(out, name)+".gcc"
            
            print name
            
            #sanefalcon            
            Parallel(n_jobs=16)(delayed(getReadStarts)(bam, id, readstarts, retropy) for id in range(1,23))
                
            #wisecondor
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


def getReadStarts(bam, id, output, retro):
	
	cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	cmd_forward += "'{print $4}' > "
	cmd_forward += "{output}.{id}.start.fwd".format(output=output, id=id)
	
	cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	cmd_reverse += "'{print ($4 + length($10) - 1)}' "
	cmd_reverse += "> {output}.{id}.start.rev".format(output=output, id=id)
	
	jobLauncher(cmd_forward)
	jobLauncher(cmd_reverse)
    
def jobLauncher(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode == 0:
		pass
        # print(stdout)
	else:
		raise Exception(stderr)

if __name__ == '__main__':
    main()