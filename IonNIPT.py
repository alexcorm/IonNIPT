#!/usr/bin/env python
from ion.plugin import * 
import os
#import sys
#import json 
import subprocess 
from joblib import Parallel, delayed
#import glob 
#import pickle
from django.template import Context, Template
from django.conf import settings

class IonNIPT(IonPlugin):
	""" IonNIPT"""
	version = "1.0"
	allow_autorun = False
	author = "cormieralexandre@gmail.com"
	envDict = dict(os.environ)
	
	def launch(self, data=None):
		print("Launch started...")
		
		# ================ GET GLOBAL PATH
		self.outputDir 		= os.environ["RESULTS_DIR"];  # The wisecondor results directory
		self.analysisDir 	= os.environ["ANALYSIS_DIR"];
		self.pluginDir		= os.environ["PLUGIN_PATH"];
		self.urlRoot 		= os.environ["URL_ROOT"]   # /output/Home/X/
		self.urlPlugin 		= os.environ["TSP_URLPATH_PLUGIN_DIR"] # /output/Home/X/plugin_out/IonWisecondor
		self.date               = os.environ["TSP_ANALYSIS_DATE"]
		
		# ================ GET INSTANCE PARAMETERS AND STORE THEM IN A LIST 
		fileCount = int(os.environ["PLUGINCONFIG__COUNT"])
		files = []
		for i in range(fileCount):
			item       	= {}
			key        	= "PLUGINCONFIG__ITEMS__"+str(i)
			barcode    	= os.environ[key+"__BARCODE"]
			sample	   	= os.environ[key+"__SAMPLE"]
			input 		= self.analysisDir +"/" + barcode + "_rawlib.bam"
			
			sample = sample.replace(' ', '_')
			
			item["sample"] 	= sample
			item["barcode"] = barcode
			item["input"] 	= input
			item["pickle"] 	= self.urlPlugin + "/" + sample + "_" + self.date +".pickle" 
			item["gcc"] 	= self.urlPlugin + "/" + sample + "_" + self.date +".gcc" 
			item["tested"] 	= self.urlPlugin + "/" + sample + "_" + self.date +".tested" 
			item["pdf"] 	= self.urlPlugin + "/" + sample + "_" + self.date +".pdf"
			
			files.append(item)
			
		# ================ LOOP ON EACH FILES AND START COMPUTATION 
		for item in files:
			# pickle
			cmd_pickle = "samtools view {bam} -q 1 | python {pluginDir}/wisecondor/consam.py -outfile {outputDir}/{sample}_{date}.pickle".format(bam=item["input"], pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			self.jobLauncher(cmd_pickle)
			
			# gcc
			cmd_gcc = "python {pluginDir}/wisecondor/gcc.py  {outputDir}/{sample}_{date}.pickle {pluginDir}/data/hg19.gccount {outputDir}/{sample}_{date}.gcc".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			self.jobLauncher(cmd_gcc)
			
			# tested
			cmd_tested = "python {pluginDir}/wisecondor/test.py {outputDir}/{sample}_{date}.gcc {pluginDir}/data/reftable {outputDir}/{sample}_{date}.tested".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			self.jobLauncher(cmd_tested)
			
			# pdf
			cmd_pdf = "python {pluginDir}/wisecondor/plot.py {outputDir}/{sample}_{date}.tested  {outputDir}/{sample}_{date}".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)	
			self.jobLauncher(cmd_pdf)
			
			# get score
			filePath = os.environ["RESULTS_DIR"] + "/" + item["sample"] + "_" + self.date + ".tested"
			item["s21"] = self.scoreOf(filePath, "21")
			item["s18"] = self.scoreOf(filePath, "18")
			item["s13"] = self.scoreOf(filePath, "13")
			
			### Saneflacon
			# prep folder
			# os.makedirs(os.path.join(self.outputDir, item["input"].split('.')[0]))
			# 
			# # readstarts
			# Parallel(n_jobs=16)(delayed(getReadStarts)(item["input"], id) for id in range(1,23))
			
			# # nucProfile
			# Parallel(n_jobs=16)(delayed(getProfiles)(item["input"].split('.')[0]), self.outputDir) for sample in os.listdir(os.path.join(self.outputDir, item["sample"])))
			# 
			# # test FF
			# cmd_ff = "predict.sh {model} nucProfile/{sample}"
			# self.jobLauncher(cmd_ff)
		
		# ================ GENERATE RESULTS HTML FROM DJANGO TEMPLATE SYSTEM
		settings.configure()
		source = open(os.environ["RUNINFO__PLUGIN__PATH"] + "/block_template.html", "r").read()
		t = Template(source)
		# Pass files arguments to the template 
		c = Context({'files': files})
		html = t.render(c)
		# Output html render 
		f = open(self.outputDir+"/resultat_block.html","w")
		f.write(html)
		f.close()

	def scoreOf(self, testedFile, chrom):
		with open(testedFile) as file:
			data 	= pickle.loads(file.read())
			zScores = data["zSmoothDict"][chrom]
			score   = -1
			try:
				score = sum(zScores) / len(zScores)
			except :
				score = -1
		return round(score,2)
	
	def jobLauncher(self, cmd):
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		stdout, stderr = p.communicate()
		if p.returncode == 0:
			print(stdout)
		else:
			raise Exception(stderr)
	
	# def getReadStarts(self, bam, id):
	# 	
	# 	retro = self.pluginDir + '/sanefalcon/retro.py'
	# 	name, extension = bam.split('.')
	# 	
	# 	cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	# 	cmd_forward += "'{print $4}' > "
	# 	cmd_forward += "{name}/{name}.{id}.start.fwd".format(name=name, id=id)
	# 	
	# 	cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	# 	cmd_reverse += "'{print ($4 + length($10) - 1)}' "
	# 	cmd_reverse += "> {name}/{name}.{id}.start.rev".format(name=name, id=id)
	# 	
	# 	process = subprocess.Popen(cmd_forward, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 	process.communicate()
	# 	process = subprocess.Popen(cmd_reverse, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 	process.communicate()
	# 	
	# def getProfiles(self, sample, dataDir):
	# 
	# 	nuclTrak = '/home/acormier/working_directory/DPNI/data/roy-nucleosome-track'
	# 	getProfile = self.pluginDir + '/getProfile.py'
	# 	
	# 	name, chr, type, strand = sample.split('.')
	# 	
	# 	if strand == 'fwd':
	# 		cmd_0 = "python {getProfile} {nuclTrak}/nucl_exR.{chr} {name}/{name}.{chr}.start.fwd 0 {name}/{name}.{chr}.fwd".format(getProfile=getProfile, nuclTrak=nuclTrak, input=input, chr=chr, sample=sample, name=name, output=output)
	# 		cmd_1 = "python {getProfile} {nuclTrak}/nucl_exR.{chr} {name}/{name}.{chr}.start.fwd 1 {name}/{name}.{chr}.ifwd".format(getProfile=getProfile, nuclTrak=nuclTrak, input=input, chr=chr, sample=sample, name=name, output=output)
	# 		
	# 		process = subprocess.Popen(cmd_0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 		process.communicate()
	# 		
	# 		process = subprocess.Popen(cmd_1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 		process.communicate()
	# 		
	# 	elif strand == 'rev':
	# 	       cmd_0 = "python {getProfile} {nuclTrak}/nucl_exR.{chr} {name}/{name}.{chr}.start.rev 1 {name}/{name}.{chr}.rev".format(getProfile=getProfile, nuclTrak=nuclTrak, input=input, chr=chr, sample=sample, name=name, output=output)
	# 	       cmd_1 = "python {getProfile} {nuclTrak}/nucl_exR.{chr} {name}/{name}.{chr}.start.rev 0 {name}/{name}.{chr}.irev".format(getProfile=getProfile, nuclTrak=nuclTrak, input=input, chr=chr, sample=sample, name=name, output=output)
	# 	       
	# 	       process = subprocess.Popen(cmd_0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 	       process.communicate()
	# 	       
	# 	       process = subprocess.Popen(cmd_1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# 	       process.communicate()
	
if __name__ == "__main__":
  PluginCLI()
