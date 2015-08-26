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

STRAND_FORWARD = True
STRAND_REVERSE = False

LEFT = True
REVERSE = False

class Fusion:
	def __init__(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence,arg_left_strand,arg_right_strand,arg_dataset_name):
		self.annotated_genes_left = None
		self.annotated_genes_right = None
		
		self.left_strand = None
		self.right_strand = None
		
		self.tested_datasets = {arg_dataset_name:True}
		self.matched_datasets = {arg_dataset_name:True}
		
		self.locations = []
		self.matches = set([self])## initial (non merged) objects used for matching
		
		self.dataset_name = arg_dataset_name
		
		self.set( \
			self.cleanup_chr_name(arg_left_chr), \
			self.cleanup_chr_name(arg_right_chr), \
			arg_left_pos, \
			arg_right_pos, \
			arg_sequence, \
			arg_transition_sequence, \
			self.find_strand_type(arg_left_strand), \
			self.find_strand_type(arg_right_strand) \
		)
	
	def get_dataset_statistics(self):
		matches = 0
		unmatches = 0
		
		for match in self.matched_datasets:
			if(match in self.tested_datasets):
				matches += 1
			else:
				unmatches += 1
		
		return (matches,unmatches)
	
	def set(self,arg_left_chr,arg_right_chr,arg_left_pos,arg_right_pos,arg_sequence,arg_transition_sequence,arg_left_strand,arg_right_strand):
		self.left_chr_str = arg_left_chr
		self.right_chr_str = arg_right_chr
		
		self.left_break_position = int(str(arg_left_pos).replace(",",""))
		self.right_break_position = int(str(arg_right_pos).replace(",",""))
		
		self.sequence = arg_sequence
		self.transition_sequence = arg_transition_sequence
		
		self.left_strand = arg_left_strand
		self.right_strand = arg_right_strand
		
		if (self.get_left_chromosome(False) > self.get_right_chromosome(False)) or ((self.get_left_chromosome(False) == self.get_right_chromosome(False)) and (self.get_left_break_position() > self.get_right_break_position())):
			self.swap_positions()
	
	def find_strand_type(self,strand_type):
		if(strand_type in [STRAND_FORWARD, STRAND_REVERSE]):
			return strand_type
		
		if isinstance(strand_type, basestring):
			strand_type = strand_type.lower()
			if strand_type in ["f","+","forward","forwards","positive","5' -> 3'"]:
				return STRAND_FORWARD
			elif strand_type in ["b","-","r","backward","backwards","reverse","negative","3' -> 5'"]:
				return STRAND_REVERSE
		
		return None
	
	def cleanup_chr_name(self,chr_name):
		"""Given the large number of fusion genes, we remove all 'chr'
		prefixes because they add 6 bytes per fusion gene. They can be
		added again using the getter functions.
		"""
		
		chr_name = chr_name.strip()
		return chr_name[3:] if chr_name[0:3] == "chr" else chr_name
	
	def add_location(self,location):
		"""Multiple locations are stored as a list. This is used in particular for merging matched fusions.
		"""
		self.locations.append(location)
	
	def get_left_position(self,indexed_chromosome=False):
		return [self.get_left_chromosome(indexed_chromosome),self.get_left_break_position()]
	
	def get_right_position(self,indexed_chromosome=False):
		return [self.get_right_chromosome(indexed_chromosome),self.get_right_break_position()]
	
	def get_left_chromosome(self,with_prefix=False):
		if(with_prefix):
			return 'chr'+self.left_chr_str
		else:
			return self.left_chr_str
	
	def get_right_chromosome(self,with_suffix=False):
		if(with_suffix):
			return 'chr'+self.right_chr_str
		else:
			return self.right_chr_str
	
	def get_left_break_position(self):
		return self.left_break_position
	
	def get_right_break_position(self):
		return self.right_break_position
	
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
		self.set( \
			self.right_chr_str, \
			self.left_chr_str, \
			self.get_right_break_position(), \
			self.get_left_break_position(), \
			self.sequence, \
			self.get_transition_sequence(), \
			self.right_strand, \
			self.left_strand \
		)
	
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
	
	def show_me(self):
		print self.__str__()
	
	def __str__(self):
		pos_left = self.get_left_position(True)
		pos_right = self.get_right_position(True)
		
		out = "Fusion (from dataset '"+self.dataset_name+"'): "+self.get_left_chromosome(True)+":"+str(pos_left[1])+"-"+self.get_right_chromosome(True)+":"+str(pos_right[1]) + "\n"
		if(self.get_annotated_genes_left()):
			out += " - annotated genes left:  "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_left()])
			if(self.left_strand == STRAND_FORWARD):
				out += " (+)"
			elif(self.left_strand == STRAND_REVERSE):
				out += " (-)"
			
			out += "\n"
		if(self.get_annotated_genes_right()):
			out += " - annotated genes right: "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_right()])
			if(self.right_strand == STRAND_FORWARD):
				out += " (+)"
			elif(self.right_strand == STRAND_REVERSE):
				out += " (-)"
			
			out += "\n"
		return out
