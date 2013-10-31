#!/usr/bin/env python



class HighThroughputFusionDetectionExperiment:
	def __init__(self,name,arg_type):
		self.name = name
		
		self.fusions_left_keys = {}
		self.fusions_right_keys = {}
		
		self.fusions_indexed_left = False
		self.fusions_indexed_right = False
		
		self.set_type(arg_type)
	
	def set_type(self,arg_type):
		if(arg_type in ["RNA","DNA"]):
			self.type = arg_type
		else:
			raise TypeError,"Incorrect type: "+str(arg_type)
	
	def get_type(self):
		return self.type
	
	def add_fusion(self,fusion):
		# Add left location
		left_chr = fusion.get_left_chromosome(False)
		
		if(not self.fusions_left_keys.has_key(left_chr)):
			self.fusions_left_keys[left_chr] = {}
		
		left_pos = fusion.get_left_break_position()
		
		if(not self.fusions_left_keys[left_chr].has_key(left_pos)):
			self.fusions_left_keys[left_chr][left_pos] = []
		self.fusions_left_keys[left_chr][left_pos].append(fusion)
		
		
		
		# Add right location
		right_chr = fusion.get_right_chromosome(False)
		
		if(not self.fusions_right_keys.has_key(right_chr)):
			self.fusions_right_keys[right_chr] = {}
		
		right_pos = fusion.get_right_break_position()
		
		if(not self.fusions_right_keys[right_chr].has_key(right_pos)):
			self.fusions_right_keys[right_chr][right_pos] = []
		self.fusions_right_keys[right_chr][right_pos].append(fusion)
	
	def count_fusions(self):
		return len(self.get_fusions())
	
	def index_fusions_left(self):
		if(not self.converted_to_genes_left()):
			print "   - Indexing (left): "+self.name
			self.fusions_indexed_left = []
			
			for left_chr in sorted(self.fusions_left_keys.keys()):
				fusions = []
				#for location in sorted(self.fusions_left_keys[left_chr]):
				#	fusions.append(self.fusions_left_keys[left_chr][location])
				for location in sorted(self.fusions_left_keys[left_chr]):
					for fusion in self.fusions_left_keys[left_chr][location]:
						fusions.append(fusion)
				
				self.fusions_indexed_left.append({"name":left_chr,"fusions":fusions})
	
	def index_fusions_right(self):
		if(not self.converted_to_genes_right()):
			print "   - Indexing (right): "+self.name
			self.fusions_indexed_right = []
			
			for right_chr in sorted(self.fusions_right_keys.keys()):
				fusions = []
				for location in sorted(self.fusions_right_keys[right_chr]):
					fusions.append(self.fusions_right_keys[right_chr][location])
				
				self.fusions_indexed_right.append({"name":right_chr,"fusions":fusions})
	
	def find_fusions(self,fusions):
		print "Number of fusions:"
		print len(fusions)
	
	def get_fusions(self):# Make iterator object?
		"""
		Returns a SORTED list of all fusions
		"""
		if(self.fusions_indexed_left == False):
			self.index_fusions_left()
		
		return self.fusions_indexed_left
	
	def get_fusions_indexed_left(self):
		if(self.fusions_indexed_left == False):
			self.index_fusions_left()
		
		return self.fusions_indexed_left
	
	def get_fusions_indexed_right(self):
		if(self.fusions_indexed_right == False):
			self.index_fusions_right()
		
		return self.fusions_indexed_right
	
	def converted_to_genes_left(self):
		return (self.fusions_indexed_left == type([]))
	
	def converted_to_genes_right(self):
		return (self.fusions_indexed_right == type([]))
	
	#def set_converted_to_genes_left(self,state):
	#	self.is_converted_to_genes = bool(state)
	
	def show_me(self):
		fusions = self.get_fusions_indexed_left()
		for _chr in fusions:
			for position in _chr["fusions"]:
				for fusion in position:
					fusion.show_me()
		print "---------------------"
	