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


import logging,fuma

from fuma import Fusion
from fuma.Fusion import AD_DIRECTION_REVERSE
from fuma.Fusion import AD_DIRECTION_FORWARD

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE


class MergedFusion:
	logger = logging.getLogger("FuMa::MergedFusion")
	
	def __init__(self):
		self.fusions = set()
	
	def add_fusion(self,arg_fusion):
		if not isinstance(arg_fusion, Fusion):
			raise Exception("MergedFusion objects can only be expanded with Fusion objects")
		else:
			len_a = len(self.fusions)
			self.fusions.update(fusion)
			if len(self.fusions) == len_a:
				raise Exception("MergedFusion is updated with one that it already contains")
	
	def locations(self):
		out = []
		for fusion in self.fusions:
			for location in fusion.locations()
				out.append(location)
		return out
	
	def get_dataset_statistics(self):
		matches = 0
		unmatches = 0
		
		for match in self.matched_datasets:
			if(match in self.tested_datasets):
				matches += 1
			else:
				unmatches += 1
		
		return (matches,unmatches)
	
	def find_strand_type(self,strand_type):
		if(strand_type in [STRAND_FORWARD, STRAND_REVERSE]):
			return strand_type
		
		if isinstance(strand_type, basestring):
			strand_type = strand_type.lower()
			if strand_type in ["f","+","forward","forwards","positive","5' -> 3'"]:
				return STRAND_FORWARD
			elif strand_type in ["b","-","r","backward","backwards","reverse","negative","3' -> 5'"]:
				return STRAND_REVERSE
			else:
				raise Exception("Unknown fusion strand: '"+strand_type+"'")
		
		return None
	
	def get_left_position(self,chromosome_with_prefix=False):
		return [self.get_left_chromosome(chromosome_with_prefix),self.get_left_break_position()]
	
	def get_right_position(self,chromosome_with_prefix=False):
		return [self.get_right_chromosome(chromosome_with_prefix),self.get_right_break_position()]
	
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
	
	def get_distance(self):
		if(self.is_interchromosomal()):
			return -1
		else:
			return self.get_right_break_position() - self.get_left_break_position()
	
	def is_interchromosomal(self):
		return (self.get_left_chromosome() != self.get_right_chromosome())
	
	def is_intrachromosomal(self):
		return (self.get_left_chromosome() == self.get_left_chromosome())
	
	def spans_a_large_gene(self):
		for gene in self.annotated_genes_left:
			if gene.is_long_gene:
				return True
		
		for gene in self.annotated_genes_right:
			if gene.is_long_gene:
				return True
		
		return False
	
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
		out = "--- MergedFusion container of size "+str(len(self.fusions))+" ---\n"
		for fusion in self.fusions:
			out += "\n"+fusion.__str__()+"\n"
		return out+"--- ---\n"
