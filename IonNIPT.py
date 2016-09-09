#!/usr/bin/env python
from ion.plugin import * 
import os
import subprocess 
from joblib import Parallel, delayed
import pickle
from django.template import Context, Template
from django.conf import settings

class IonNIPT(IonPlugin):
	""" IonNIPT"""
	version = "1.0"
	allow_autorun = False
	author = "cormieralexandre@gmail.com"
	envDict = dict(os.environ)
	
	def launch(self, data=None):
		print "Launch started..."
		
		# ================ GET GLOBAL PATH
		self.outputDir 		= os.environ["RESULTS_DIR"];  # The plugin results directory
		self.analysisDir 	= os.environ["ANALYSIS_DIR"];
		self.pluginDir		= os.environ["PLUGIN_PATH"];
		self.urlRoot 		= os.environ["URL_ROOT"]   # /output/Home/X/
		self.urlPlugin 		= os.environ["TSP_URLPATH_PLUGIN_DIR"] # /output/Home/X/plugin_out/IonWisecondor
		self.date               = os.environ["TSP_ANALYSIS_DATE"]
		
		self.retro = os.path.join(self.pluginDir, 'sanefalcon/retro.py')
		self.getProfile = os.path.join(self.pluginDir, 'sanefalcon/getProfile.py')
		self.predict = os.path.join(self.pluginDir, 'sanefalcon/predict.sh')
		self.defrag = os.path.join(self.pluginDir, 'wisecondor/defrag.py')
		self.nuclTrack = os.path.join(self.pluginDir, 'data')
		self.trainModel = os.path.join(self.pluginDir, 'data/trainModel-brest.model')
		self.male = os.path.join(self.pluginDir, 'data/male-defrag')
		self.female = os.path.join(self.pluginDir, 'data/female-defrag')
		self.sf = 0.688334125062
		self.percYonMales = 0.00146939199267
		
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
			print "Analyse {sample}".format(sample=item["sample"])
			
			# ======= Wisecondor
			
			print "*** Detect trisomy using Wisecondor ***"
			print "Get pickle"
			# pickle
			cmd_pickle = "samtools view {bam} -q 1 | python {pluginDir}/wisecondor/consam.py -outfile {outputDir}/{sample}_{date}.pickle".format(bam=item["input"], pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			jobLauncher(cmd_pickle)
			
			print "Make GC correction"
			# gcc
			cmd_gcc = "python {pluginDir}/wisecondor/gcc.py  {outputDir}/{sample}_{date}.pickle {pluginDir}/data/hg19.gccount {outputDir}/{sample}_{date}.gcc".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			jobLauncher(cmd_gcc)
			
			print "Test"
			# tested
			cmd_tested = "python {pluginDir}/wisecondor/test.py {outputDir}/{sample}_{date}.gcc {pluginDir}/data/reftable {outputDir}/{sample}_{date}.tested".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)
			jobLauncher(cmd_tested)
			
			print "Generate the pdf file"
			# pdf
			cmd_pdf = "python {pluginDir}/wisecondor/plot.py {outputDir}/{sample}_{date}.tested  {outputDir}/{sample}_{date}".format(pluginDir = self.pluginDir, outputDir = self.outputDir, sample = item["sample"], date= self.date)	
			jobLauncher(cmd_pdf)
			
			# get score
			filePath = os.environ["RESULTS_DIR"] + "/" + item["sample"] + "_" + self.date + ".tested"
			item["s21"] = self.scoreOf(filePath, "21")
			item["s18"] = self.scoreOf(filePath, "18")
			item["s13"] = self.scoreOf(filePath, "13")
			
			# ======= Saneflacon
			
			print "*** Get foetal fraction using Sanefalcon ***"
			
			# prep folder
			os.makedirs(os.path.join(self.outputDir, item["barcode"]))
			barcode_folder = os.path.join(self.outputDir, item["barcode"])
			
			# readstarts
			print "Compute readstarts"
			Parallel(n_jobs=16)(delayed(getReadStarts)(item["input"], id, barcode_folder, self.retro) for id in range(1,23))
			
			# nucProfile
			print "Get nucleosome profile"
			Parallel(n_jobs=16)(delayed(getProfiles)(os.path.join(barcode_folder, readstarts), self.getProfile, self.nuclTrack, barcode_folder) for readstarts in os.listdir(barcode_folder))
			
			# test FF
			print "Compute the foetal fraction"
			cmd_ff = "bash {predict} {trainModel} {barcode_folder}/{barcode}".format(predict=self.predict, trainModel=self.trainModel, barcode_folder=barcode_folder, barcode=item["barcode"])
			jobLauncher(cmd_ff)
			
			ff = open(os.path.join(barcode_folder, item["barcode"])+'.ff','r')
			for line in ff:
				elem = line.split()
				if elem[0] == 'Fetal':
					item["ff_sanefalcon"] = round(float(elem[2]), 2)
			
			# ======= Defrag
			
			src_gcc = os.path.join(self.outputDir, item["sample"] + "_" + self.date +".gcc")
			src_pickle = os.path.join(self.outputDir, item["sample"] + "_" + self.date +".pickle")
			
			dest_gcc = os.path.join(barcode_folder, item["barcode"]+'.gcc')
			dest_pickle = os.path.join(barcode_folder, item["barcode"]+'.pickle')
			
			os.symlink(src_gcc, dest_gcc)
			os.symlink(src_pickle, dest_pickle)
			
			cmd_defrag = 'python {defrag} {male} {female} {test} --scalingFactor {sf} --percYonMales {percYmales} {graph}'.format(defrag=self.defrag, male=self.male, female=self.female , test=barcode_folder, sf=self.sf, percYmales=self.percYonMales, graph=item["barcode"])
			p = subprocess.Popen(cmd_defrag, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			stdout, stderr = p.communicate()
			if p.returncode == 0:
				for line in stdout.splitlines():
					if line.startswith('Ion'):
						barcode, ff, ffWholeY, gender, reads, cluster, percReadsY = line.split()
						item["ff_defrag"] = round(float(ff), 2)
						item["sex"] = gender
						item["cluster"] = cluster
			else:
				raise Exception(stderr)
			
		
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
	
def jobLauncher(cmd):
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if p.returncode == 0:
		print(stdout)
	else:
		raise Exception(stderr)
	
def getReadStarts(bam, id, output, retro): #remove from class, joblib error
	name = os.path.basename(bam).split('.')[0][:-7]
	
	cmd_forward = "samtools view {bam} chr{id} -F 20 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	cmd_forward += "'{print $4}' > "
	cmd_forward += "{output}/{name}.{id}.start.fwd".format(output=output, name=name, id=id)
	
	cmd_reverse = "samtools view {bam} chr{id} -f 16 -F 4 -q 1 | python {retro} | awk ".format(bam=bam, retro=retro, id=id)
	cmd_reverse += "'{print ($4 + length($10) - 1)}' "
	cmd_reverse += "> {output}/{name}.{id}.start.rev".format(output=output, name=name, id=id)
	
	jobLauncher(cmd_forward)
	jobLauncher(cmd_reverse)
		
def getProfiles(readstarts, getProfile, nuclTrak, barcode_folder):
	name, chr, type, strand = os.path.basename(readstarts).split('.')
	
	if strand == 'fwd':
		cmd_0 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {barcode_folder}/{name}.{chr}.start.fwd 0 {barcode_folder}/{name}.{chr}.fwd".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, barcode_folder=barcode_folder)
		cmd_1 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {barcode_folder}/{name}.{chr}.start.fwd 1 {barcode_folder}/{name}.{chr}.ifwd".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, barcode_folder=barcode_folder)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)
		
	elif strand == 'rev':
		cmd_0 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {barcode_folder}/{name}.{chr}.start.rev 1 {barcode_folder}/{name}.{chr}.rev".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, barcode_folder=barcode_folder)
		cmd_1 = "python {getProfile} {nuclTrak}/nuclTrack.{chr} {barcode_folder}/{name}.{chr}.start.rev 0 {barcode_folder}/{name}.{chr}.irev".format(getProfile=getProfile, nuclTrak=nuclTrak, chr=chr, name=name, barcode_folder=barcode_folder)
		
		jobLauncher(cmd_0)
		jobLauncher(cmd_1)
	
if __name__ == "__main__":
  PluginCLI()
