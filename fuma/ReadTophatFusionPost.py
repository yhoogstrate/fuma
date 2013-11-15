#!/usr/bin/env python



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class ReadTophatFusionPost(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		self.filename = arg_filename
		
		self.parse()
	
	def flush(self):
		self.chr_1 = False
		self.chr_2 = False
		
		self.break_1 = False
		self.break_2 = False
		
		self.seq = False
		self.insert_seq = False
	
	def parse_line_type_0(self,line):
		line = line.strip().split(" ")
		_chr = line[1].split("-")
		
		self.chr_1 = _chr[0]
		self.chr_2 = _chr[1]
		
		self.break_1 = int(line[2])
		self.break_2 = int(line[3])
		
		self.left_strand = line[4][0]
		self.right_strand = line[4][1]
	
	def parse_line_type_1(self,line):
		line = line.strip().split(" ")
		if(len(line) > 0):
			self.seq = line[0]
		
	def parse_line_type_2(self,line):
		line = line.strip().split(" ")
		if(len(line) > 1):
			self.seq += line[1]
	
	def parse(self):
		self.flush()
		
		with open(self.filename,"r") as fh:
			i = 0
			for line in fh:
				line_type = i % 6
				if(line_type == 0):
					self.parse_line_type_0(line)
				elif(line_type == 1):
					self.parse_line_type_1(line)
				elif(line_type == 2):
					self.parse_line_type_2(line)
					self.add_fusion(Fusion(self.chr_1,self.chr_2,self.break_1,self.break_2,self.seq,self.insert_seq,self.left_strand,self.right_strand,self.name))
				
				i += 1

