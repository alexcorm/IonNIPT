#!/usr/bin/env python
# -*- coding: utf-8 -*-

#================================================================#

import os
import sys
import re
import json
import subprocess

#================================================================#

def main():
    
    # Y linked genes
    
    Y = {
        "HSFY1" : "chrY:20708557-20750849",
        #"HSFY2" : "chrY:20893326-20990548",
        "BPY2" : "chrY:25119966-25151612",
        "BPY2B" : "chrY:26753707-26785354",
        "BPY2C" : "chrY:27177048-27208695",
        "XKRY " : "chrY:19880860-19889280",
        "PRY" : "chrY:24636544-24660784",
        "PRY2" : "chrY:24217903-24242154"
        }
    
    path_run = sys.argv[1] # folder where bam (and json files for Ion Torrent  server) is/are located
    
    # ### Ion Torrent method
    # # load run info
    # dpniDict = {}
    # 
    # json_file = '%s/ion_params_00.json' % path_run
    # json_load = json.load(open(json_file))
    # runname = json_load['exp_json']['log']['runname'] #.split('-DANNI')[0]
    # # runname ='-'.join(runname.split('-')[0:2])
    # 
    # dpniDict[runname] = []
    # 
    # #get sample and barcode name
    # for sample, barcode in json_load['experimentAnalysisSettings']['barcodedSamples'].items():
    #     
    #     sample_name = sample.replace(' ', '_') # just in case...
    #     barcode_name = barcode['barcodeSampleInfo'].keys()[0]
    #     
    #     dpniDict[runname].append([barcode_name, sample_name])
    # 
    # for run, design in dpniDict.items():
    #     for sample in design:
    #         
    #         barcode_id = sample[0]
    #         sample_name = sample[1]
    #         name = run + "_" + sample_name
    #         bam = os.path.join(path_run, barcode_id)+"_rawlib.bam"
    #         
    #         Yl_reads = []
    #         Y_reads = []
    #         
    #         #coverage Y linked genes
    #         for gene,coord in Y.items():
    #             cmd = "samtools view -c {bam} {coord}".format(bam = bam, coord = coord)
    #             process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #             stdout, stderr = process.communicate()
    #             if process.returncode != 0:
    #                 raise Exception(stderr)
    #             else:
    #                 Yl_reads.append(int(stdout[:-1]))
    #                 
    #         cmd = "samtools view -c {bam} {coord}".format(bam = bam, coord = "chrY")
    #         process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #         stdout, stderr = process.communicate()
    #         if process.returncode != 0:
    #             raise Exception(stderr)
    #         else:
    #             Y_reads.append(int(stdout[:-1]))
    #         
    #         print name, sum(Yl_reads), Y_reads[0], float(sum(Yl_reads))*100/float(Y_reads[0]), Yl_reads
             
    ### cluster method
    for bam in os.listdir(path_run):
        if bam.endswith(".bam") or bam.endswith(".cram"):
            name = os.path.basename(bam).split('.')[0]
            bam_path = os.path.join(path_run, bam)
            
            Yl_reads = []
            Y_reads = []
            
            for gene,coord in Y.items():
                cmd = "samtools view -c {bam} {coord}".format(bam = bam_path, coord = coord)
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                if process.returncode != 0:
                    raise Exception(stderr)
                else:
                    Yl_reads.append(int(stdout[:-1]))
                    
            cmd = "samtools view -c {bam} {coord}".format(bam = bam_path, coord = "chrY")
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                raise Exception(stderr)
            else:
                Y_reads.append(int(stdout[:-1]))
            
            if int(Y_reads[0]) != 0:
                print name, sum(Yl_reads), Y_reads[0], float(sum(Yl_reads))*100/float(Y_reads[0]), Yl_reads
            else:
                print name, sum(Yl_reads), Y_reads[0], 0, Yl_reads

if __name__ == '__main__':
    main()