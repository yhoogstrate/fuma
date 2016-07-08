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

from .Fusion import Fusion
from fuma.Fusion import AD_DIRECTION_REVERSE
from fuma.Fusion import AD_DIRECTION_FORWARD

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE


class MergedFusion:
	logger = logging.getLogger("FuMa::MergedFusion")
	
	def __init__(self):
		self.fusions = set()
		
		annotated_genes_left = None
		annotated_genes_right = None
	
	def __len__(self):
		return len(self.fusions)
	
	def add_fusion(self,arg_fusion):
		if not isinstance(arg_fusion, Fusion):
			raise Exception("MergedFusion objects can only be expanded with Fusion objects and not with: "+arg_fusion.__class__.__name__)
		else:
			len_a = len(self.fusions)
			self.fusions.add(arg_fusion)
			new_len = len(self.fusions)
			if new_len == len_a:
				raise Exception("MergedFusion is updated with one that it already contains")
			
			#if new_len == 2:
			#	self.logger.debug("Merged fusion genes")
			#elif new_len > 2:
			#	self.logger.debug("Merged fusion gene to size="+str(new_len))

	
	## After reconsidering, this function should most likely not be necesairy
	#def merge(self,arg_merged_fusion):
	#	if not isinstance(arg_merged_fusion, MergedFusion):
	#		raise Exception("MergedFusion objects can only be merged with other MergedFusion objects and not with: "+arg_fusion.__class__.__name__)
	#	
	#	len_a = len(self.fusions)
	#	
	#	for fusion_gene in arg_merged_fusion.fusions:
	#		self.fusions.add(fusion_gene)
	#	
	#	if len(self.fusions) == len_a:
	#		raise Exception("MergedFusion is merged with one that contains no new ones")
	
	def locations(self):
		out = []
		for fusion in self.fusions:
			for location in fusion.locations():
				out.append(location)
		return out
	
	def get_left_strand(self):
		strands = []
		
		for fusion in self.fusions:
			strands.append(fusion.get_left_strand())
		
		strands = list(set(strands))
		
		if len(strands) != 1:
			raise Exception("Inconsistent left-strands - if this happens either the request is unneccesary or the merge of fusions was invalid")
		else:
			return strands[0]
	
	def get_right_strand(self):
		strands = []
		
		for fusion in self.fusions:
			strands.append(fusion.get_right_strand())
		
		strands = list(set(strands))
		
		if len(strands) != 1:
			raise Exception("Inconsistent right-strands - if this happens either the request is unneccesary or the merge of fusions was invalid")
		else:
			return strands[0]
	
	def acceptor_donor_direction(self):
		directions = []
		
		for fusion in self.fusions:
			directions.append(fusion.acceptor_donor_direction())
		
		directions = list(set(directions))
		
		if len(directions) != 1:
			raise Exception("Inconsistent acceptor donor directions - if this happens either the request is unneccesary or the merge of fusions was invalid")
		else:
			return directions[0]
	
	def spans_a_large_gene(self):
		for gene in self.annotated_genes_left:
			if gene.is_long_gene:
				return True
		
		for gene in self.annotated_genes_right:
			if gene.is_long_gene:
				return True
		
		return False
	
	def get_annotated_genes_left2(self):
		if(not self.has_annotated_genes()):
			raise Exception("Requested empty gene list")
		else:
			return self.annotated_genes_left
	
	def get_annotated_genes_right2(self):
		if(not self.has_annotated_genes()):
			raise Exception("Requested empty gene list")
		else:
			return self.annotated_genes_right
	
	def has_annotated_genes(self):
		for fusion in self.fusions:
			if not fusion.has_annotated_genes():
				return False
		return True
	
	def show_me(self):
		print(self.__str__())
	
	def __str__(self):
		out = "--- MergedFusion container of size "+str(len(self.fusions))+" ---"
		n = len(out)
		for fusion in self.fusions:
			out += "\n"+fusion.__str__()+"\n"
		return out+("-"*n)+"\n"
