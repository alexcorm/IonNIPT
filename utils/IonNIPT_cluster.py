#!/usr/bin/env python

import argparse
import os
import subprocess 
from joblib import Parallel, delayed
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def getArgs():
	parser = argparse.ArgumentParser(description="",version="1.0.0")
	parser.add_argument('-b',dest="bamDir",type=str,required=True,help='BAM folder')
	parser.add_argument('-p',dest="pluginDir",type=str,required=True,help='Path to IonNIPT folder')
	parser.add_argument('-o',dest="outDir",type=str,required=True,help='Output folder')
	arg = parser.parse_args()
	
	return arg

def main(args):

	### get BAM files (allows parallelisation)
	print('Get BAM list...')
	bam_lst = []
	for bam in os.listdir(args.bamDir):
		if bam.endswith(".bam"):
			bam_lst.append(bam)
	print('Done')
	
	### Make Wisecondor analysis (in //)
	print('WiseCondor analysis running...')
	Parallel(n_jobs=16)(delayed(wisecondor)(bam, args.bamDir, args.outDir, args.pluginDir) for bam in bam_lst)
	print('Done')
	
	### Make Sanefalcon analysis
	
	print('Start Sanefalcon analysis...')
	readstarts = os.path.join(args.outDir, 'readstarts')
	os.makedirs(readstarts)
	
	## Get readstarts for each bam, in // per chromosome
	print('First step, get readstarts')
	for bam in bam_lst:	
		Parallel(n_jobs=16)(delayed(getReadStarts)(bam, id, args.bamDir, readstarts, args.pluginDir) for id in range(1,23))
	print('Readstarts: OK')
	
	## Get nucleosome profiles
	profiles = os.path.join(args.outDir, 'profiles')
	os.makedirs(profiles)
	
	print('Second step, get nucleosomes profiles')
	Parallel(n_jobs=16)(delayed(getProfiles)(readstart, args.pluginDir, readstarts, profiles) for readstart in os.listdir(readstarts))
	print('Profiles: OK')	
	## Get foetal fraction
	
	predict = os.path.join(args.pluginDir, 'sanefalcon/predict.sh')
	trainModel = os.path.join(args.pluginDir, 'data/trainModel_BorBreCoc_defrag-rassf1.model')
	print('Final step, get ff with Sanefalcon')
	ff_sanefalcon = {}
	for bam in bam_lst:
		name = os.path.basename(bam).split('.')[0]
		bam = os.path.join(profiles, name)
		
		ff_cmd = "bash {predict} {trainModel} {bam}".format(predict=predict, trainModel=trainModel, bam=bam)
		jobLauncher(ff_cmd)
		
		ff = open(os.path.join(bam + '.ff'),'r')
		for line in ff:
			elem = line.split()
			if elem[0] == 'Fetal':
				ff_sanefalcon[name] = round(float(elem[2]), 2)
	
	for key, value in sorted(ff_sanefalcon.items()):
		print key, value
		
	print('Done')
	
	#### defrag
	
	defrag = os.path.join(args.pluginDir, 'wisecondor/defrag.py')
	male = os.path.join(args.pluginDir, 'data/male-defrag')
	female = os.path.join(args.pluginDir, 'data/female-defrag')
	sf = 0.688334125062
	percYonMales = 0.00146939199267
	
	###
	
	print('Defrag analysis...')
	defrag_cmd = 'python {defrag} {male} {female} {test} --scalingFactor {sf} --percYonMales {percYmales} {graph} > defrag.res'.format(defrag=defrag, male=male, female=female , test=args.outDir, sf=sf, percYmales=percYonMales, graph="defrag")
	p = subprocess.Popen(defrag_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		raise Exception(stderr)
	
	print('Done')

def wisecondor(bam, input, output, plugin):
	
	### exe + data
	consampy = os.path.join(plugin, 'wisecondor/consam.py')
	testpy = os.path.join(plugin, 'wisecondor/test.py')
	plotpy = os.path.join(plugin, 'wisecondor/plot.py')
	gccpy = os.path.join(plugin, 'wisecondor/gcc.py')
	gccount = os.path.join(plugin, 'data/hg19.gccount')
	reftable =  os.path.join(plugin, 'data/reftable')
	
	### path + output names
	name = os.path.basename(bam).split('.')[0]
	bam = os.path.join(input, bam) # -> crush bam into new bam var, full path
	
	pickle = os.path.join(output, name + '.pickle')
	gcc = os.path.join(output, name + '.gcc')
	tested = os.path.join(output, name + '.tested')
	pdf = os.path.join(output, name)
		
	### wisecondor analysis
	pickle_cmd = "samtools view {bam} -q 1 | python {consampy} -outfile {pickle}".format(bam = bam, consampy = consampy, pickle = pickle)
	jobLauncher(pickle_cmd)
	
	gcc_cmd = "python {gccpy} {pickle} {gccount} {gcc}".format(gccpy = gccpy, pickle = pickle, gccount = gccount, gcc = gcc)
	jobLauncher(gcc_cmd)
	
	tested_cmd = "python {testpy} {gcc} {reftable} {tested}".format(testpy =  testpy, gcc = gcc, reftable = reftable, tested = tested)
	jobLauncher(tested_cmd)
	
	pdf_cmd = "python {plotpy} {tested} {pdf}".format(plotpy = plotpy, tested = tested, pdf = pdf)
	jobLauncher(pdf_cmd)

def getReadStarts(bam, id, input, output, plugin):
	
	name = os.path.basename(bam).split('.')[0]
	bam = os.path.join(input, bam)
	
	forward = os.path.join(output, name + '.' + str(id) + '.start.fwd')
	reverse = os.path.join(output, name + '.' + str(id) + '.start.rev')
	
	###
	
	retro = os.path.join(plugin, 'sanefalcon/retro.py')
	
	###
	
	cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam = bam, retro = retro, id = id)
	cmd_forward += "'{print $4}' > "
	cmd_forward += "{forward}".format(forward = forward)
	
	cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam = bam, retro = retro, id = id)
	cmd_reverse += "'{print ($4 + length($10) - 1)}' "
	cmd_reverse += "> {reverse}".format(reverse = reverse)
	
	jobLauncher(cmd_forward)
	jobLauncher(cmd_reverse)

def getProfiles(readstart, plugin, input, output):
	name, chr, type, strand = os.path.basename(readstart).split('.')
	
	###
	
	getProfile = os.path.join(plugin, 'sanefalcon/getProfile.py')
	nuclTrack = os.path.join(plugin, 'data')
	
	###
	
	if strand == 'fwd':
		cmd_0 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.fwd 0 {output}/{name}.{chr}.fwd".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		cmd_1 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.fwd 1 {output}/{name}.{chr}.ifwd".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)
		
	elif strand == 'rev':
		cmd_0 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.rev 1 {output}/{name}.{chr}.rev".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		cmd_1 = "python {getProfile} {nuclTrack}/nuclTrack.{chr} {input}/{name}.{chr}.start.rev 0 {output}/{name}.{chr}.irev".format(getProfile=getProfile, nuclTrack=nuclTrack, chr=chr, name=name, input=input, output=output)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)

def jobLauncher(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		raise Exception(stderr)

if __name__ == '__main__':
	args = getArgs()
	main(args)