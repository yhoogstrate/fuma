#!/usr/bin/env python



class Fusion:
	def __init__(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence):
		self.annotated_genes_left = False
		self.annotated_genes_right = False
		self.set(arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence)
	
	def set(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence):
		self.left_chr_str = arg_left_chr
		self.right_chr_str = arg_right_chr
		
		left_key = arg_left_chr.replace("chr","")
		right_key = arg_right_chr.replace("chr","")
		
		if(left_key == self.left_chr_str):
			self.left_chr_str = "chr"+self.left_chr_str
		if(right_key == self.right_chr_str):
			self.right_chr_str = "chr"+self.right_chr_str
		
		if(left_key.isdigit()):
			self.left_chr_key = int(left_key)
		else:
			self.left_chr_key = left_key
		
		
		if(right_key.isdigit()):
			self.right_chr_key = int(right_key)
		else:
			self.right_chr_key = right_key
		
		self.left_break_position = int(arg_left_pos)
		self.right_break_position = int(arg_right_pos)
		
		self.sequence = arg_sequence
		self.transition_sequence = arg_transition_sequence
		
		if((self.get_left_chromosome(False) > self.get_right_chromosome(False)) or ((self.get_left_chromosome(False) == self.get_right_chromosome(False)) and (self.get_left_break_position() > self.get_right_break_position()))):
			self.swap_positions()
	
	def get_left_position(self,indexed_chromosome=False):
		return [self.get_left_chromosome(indexed_chromosome),self.get_left_break_position()]
	
	def get_right_position(self,indexed_chromosome=False):
		return [self.get_right_chromosome(indexed_chromosome),self.get_right_break_position()]
	
	def get_left_chromosome(self,with_suffix=True):
		if(with_suffix):
			return self.left_chr_str
		else:
			return self.left_chr_key
	
	def get_right_chromosome(self,with_suffix=True):
		if(with_suffix):
			return self.right_chr_str
		else:
			return self.right_chr_key
	
	def get_left_break_position(self):
		return self.left_break_position
	
	def get_right_break_position(self):
		return self.right_break_position
	
	def get_sequence(self):
		return self.sequence
	
	def get_transition_sequence(self):
		return self.transition_sequence
	
	def swap_positions(self):
		self.set(self.get_right_chromosome(),self.get_left_chromosome(),self.get_right_break_position(),self.get_left_break_position(),self.get_sequence(),self.get_transition_sequence())
	
	def annotate_genes_left(self,gene_names):
		self.annotated_genes_left = gene_names
		
	def annotate_genes_right(self,gene_names):
		self.annotated_genes_right = gene_names
	
	def get_annotated_genes_left(self):
		return self.annotated_genes_left
	
	def get_annotated_genes_right(self):
		return self.annotated_genes_right
	
	def show_me(self):
		print "Fusion:"
		print " - left:",self.get_left_position(),"\n - right:",self.get_right_position(),"\n - annotated genes left:",self.get_annotated_genes_left(),"\n - annotated genes right:",self.get_annotated_genes_right()

