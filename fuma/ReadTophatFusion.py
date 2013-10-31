#!/usr/bin/env python



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class ReadTophatFusion(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		if(len(line) > 0):
			q = line
			line = line.split("@")
			line[0] = line[0].split("\t")
			line[2] = line[2].strip().split(" ")
			line[3] = line[3].strip().split(" ")
			
			chromosomes = line[0][0].split("-")
			
			if(len(line[3]) > 1):
				sequence = line[2][0] + line[3][1]
			else:
				sequence = False
			
			self.add_fusion(Fusion(chromosomes[0],chromosomes[1],line[0][1],line[0][2],sequence,False))
	
	def parse(self):
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
