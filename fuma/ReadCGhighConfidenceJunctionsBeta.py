#!/usr/bin/env python



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class ReadCGhighConfidenceJunctionsBeta(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"DNA")
		
		self.filename = arg_filename
		
		self.parse_left_chr_column = -1
		self.parse_right_chr_column = -1
		
		self.parse_left_pos_column = -1
		self.parse_right_pos_column = -1
		
		self.parse_sequence_column = -1
		self.parse_transition_sequence_column = -1
		
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		if(len(line) > 0):
			if(line[0] != "#"):
				if(line[0] == ">"):
					self.parse_line__header(line)
				else:
					self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line[1:]
		line = line.split("\t")
		
		self.parse_left_chr_column = line.index("LeftChr")
		self.parse_right_chr_column = line.index("RightChr")
		
		self.parse_left_pos_column = line.index("LeftPosition")
		self.parse_right_pos_column = line.index("RightPosition")
		
		self.parse_sequence_column = line.index("AssembledSequence")
		self.parse_transition_sequence_column = line.index("TransitionSequence")
		
		self.parse_left_strand = line.index("LeftStrand")
		self.parse_right_strand = line.index("RightStrand")
	
	def parse_line__fusion(self,line):
		line = line.split("\t")
		
		self.add_fusion(Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],line[self.parse_left_pos_column],line[self.parse_right_pos_column],line[self.parse_sequence_column],line[self.parse_transition_sequence_column],line[self.parse_left_strand],line[self.parse_right_strand],self.name))
	
	def parse(self):
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
