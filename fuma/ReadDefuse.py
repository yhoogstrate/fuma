#!/usr/bin/env python



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class ReadDefuse(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		if(self.parse_header):
			self.parse_line__header(line)
		else:
			self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line.strip().split("\t")
		
		self.parse_left_chr_column = line.index("gene_chromosome1")
		self.parse_right_chr_column = line.index("gene_chromosome2")
		
		self.parse_left_pos_column = line.index("genomic_break_pos1")
		self.parse_right_pos_column = line.index("genomic_break_pos2")
		
		self.parse_sequence_column = line.index("splitr_sequence")
		
		self.parse_header = False
	
	def parse_line__fusion(self,line):
		line = line.strip().split("\t")
		self.add_fusion(Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],line[self.parse_left_pos_column],line[self.parse_right_pos_column],line[self.parse_sequence_column],False))
	
	def parse(self):
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
