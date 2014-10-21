#!/usr/bin/env python

"""[License: GNU General Public License v3 (GPLv3)]
 
 This file is part of FuMa.
 
 FuMa is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 FuMa is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

 Documentation as defined by:
 <http://epydoc.sourceforge.net/manual-fields.html#fields-synonyms>
"""

class Fusion:
	def __init__(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence,arg_left_strand,arg_right_strand,arg_dataset_name):
		self.annotated_genes_left = None
		self.annotated_genes_right = None
		
		#self.left_start = False
		#self.right_start = False
		
		self.locations = []
		
		self.dataset_name = arg_dataset_name
		#self.datasets = []												#@create list with dataset ids - strings consume quite some space...
		
		self.set(arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence,arg_left_strand,arg_right_strand)
	
	def set(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence,arg_left_strand,arg_right_strand):
		self.left_chr_str = arg_left_chr
		self.right_chr_str = arg_right_chr
		
		if(arg_left_chr[0:3] == "chr"):
			 arg_left_chr = arg_left_chr.replace("chr","")
		if(arg_right_chr[0:3] == "chr"):
			arg_right_chr = arg_right_chr.replace("chr","")
		
		left_key = arg_left_chr
		right_key = arg_right_chr
		
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
		
		self.left_break_position = int(str(arg_left_pos).replace(",",""))
		self.right_break_position = int(str(arg_right_pos).replace(",",""))
		
		self.sequence = arg_sequence
		self.transition_sequence = arg_transition_sequence
		
		if(arg_left_strand in ("f","F","+","forward","forwards")):
			self.left_strand = "+"
		elif(arg_left_strand in ("b","B","-","r","R","backward","backwards","reverse")):
			self.left_strand = "-"
		else:
			raise TypeError,"unknown strand: ",arg_left_strand
		
		if(arg_right_strand in ("f","F","+","forward","forwards","positive")):
			self.right_strand = "+"
		elif(arg_right_strand in ("b","B","-","r","R","backward","backwards","reverse","negative")):
			self.right_strand = "-"
		else:
			raise TypeError,"unknown strand: ",arg_right_strand
		
		if((self.get_left_chromosome(False) > self.get_right_chromosome(False)) or ((self.get_left_chromosome(False) == self.get_right_chromosome(False)) and (self.get_left_break_position() > self.get_right_break_position()))):
			self.swap_positions()
	
	def add_location(self,location):
		"""
		For multiple locations.. tricky stuff
		"""
		self.locations.append(location)
	
	def get_left_position(self,indexed_chromosome=False):
		return [self.get_left_chromosome(indexed_chromosome),self.get_left_break_position()]
	
	def get_right_position(self,indexed_chromosome=False):
		return [self.get_right_chromosome(indexed_chromosome),self.get_right_break_position()]
	
	def get_left_chromosome(self,with_prefix=True):
		if(with_prefix):
			return self.left_chr_str
		else:
			return self.left_chr_key
	
	def get_right_chromosome(self,with_suffix=True):
		if(with_suffix):
			return self.right_chr_str
		else:
			return self.right_chr_key
	
	def get_left_strand(self):
		return self.left_strand
	
	def get_right_strand(self):
		return self.right_strand
	
	def get_left_break_position(self):
		return self.left_break_position
	
	def get_right_break_position(self):
		return self.right_break_position
	
	def get_sequence(self):
		return self.sequence
	
	def get_transition_sequence(self):
		return self.transition_sequence
	
	def get_distance(self):
		if(self.is_interchromosomal()):
			return -1
		else:
			return self.get_right_break_position() - self.get_left_break_position()
	
	def is_interchromosomal(self):
		return (self.get_left_chromosome() != self.get_right_chromosome())
	
	def is_intrachromosomal(self):
		return (self.get_left_chromosome() == self.get_left_chromosome())
	
	def swap_positions(self):
		self.set(self.get_right_chromosome(),self.get_left_chromosome(),self.get_right_break_position(),self.get_left_break_position(),self.get_sequence(),self.get_transition_sequence(),self.get_right_strand(),self.get_left_strand())
	
	def annotate_genes_left(self,gene_names):
		self.annotated_genes_left = gene_names
		
	def annotate_genes_right(self,gene_names):
		self.annotated_genes_right = gene_names
	
	def get_annotated_genes_left(self,name_indexed = False):
		if(not name_indexed):
			if(not self.annotated_genes_left):
				return []
			else:
				return self.annotated_genes_left
		else:
			index = {}
			
			for gene in self.get_annotated_genes_left():
				gene_name = str(gene)
				if(not index.has_key(gene_name)):
					index[gene_name] = []
				index[gene_name].append(gene)
			
			return index
	
	def get_annotated_genes_right(self,name_indexed = False):
		if(not name_indexed):
			if(not self.annotated_genes_right):
				return []
			else:
				return self.annotated_genes_right
		else:
			index = {}
			for gene in self.get_annotated_genes_right():
				gene_name = str(gene)
				if(not index.has_key(gene_name)):
					index[gene_name] = []
				
				index[gene_name].append(gene)
			
			return index
	
	def get_dataset_name(self):
		"""@todo get_dataset_id(self):
		"""
		return self.dataset_name
	
	def show_me(self):
		pos_left = self.get_left_position(True)
		pos_right = self.get_right_position(True)
		
		print "Fusion (from "+self.get_dataset_name()+"): "+pos_left[0]+":"+str(pos_left[1])+"-"+pos_right[0]+":"+str(pos_right[1])
		if(self.get_annotated_genes_left()):
			print " - annotated genes left:  "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_left()])
		if(self.get_annotated_genes_right()):
			print " - annotated genes right: "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_right()])
