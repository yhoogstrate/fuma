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

STRAND_FORWARD      = True
STRAND_REVERSE      = False
STRAND_UNDETERMINED = None

# Acceptor-Donor direction equal to lexicographical order?
AD_DIRECTION_FORWARD = True
AD_DIRECTION_REVERSE = False

class Fusion:
	def __init__(self, \
	   arg_left_chr, \
	   arg_right_chr, \
	   
	   arg_left_pos, \
	   arg_right_pos, \
	   
	   arg_left_strand, \
	   arg_right_strand, \
	   
	   arg_dataset_name, \
	   
	   arg_uid, \
	   arg_auto_set_acceptor_donor_direction# If a particular Fusion Gene can either be 5' -> 3' or 5' <- 3', this has to be set to False
	   ):
		self.annotated_genes_left = None
		self.annotated_genes_right = None
		
		self.left_strand = None
		self.right_strand = None
		
		self.tested_datasets = {arg_dataset_name:True}
		self.matched_datasets = {arg_dataset_name:True}
		
		self.matches = set([self])## initial (non merged) objects used for matching
		
		#@todo use pointer to original dataset?
		self.dataset_name = arg_dataset_name
		
		self.acceptor_donor_direction = None
		
		self.set( \
			self.cleanup_chr_name(arg_left_chr), \
			self.cleanup_chr_name(arg_right_chr), \
			arg_left_pos, \
			arg_right_pos, \
			self.find_strand_type(arg_left_strand), \
			self.find_strand_type(arg_right_strand), \
			arg_uid, \
			arg_auto_set_acceptor_donor_direction
		)
	
	def prepare_deletion(self):
		del(self.matches)
	
	def locations(self):
		out = []
		for match in self.matches:
			location = {
				'left':[match.get_left_chromosome(),  match.get_left_break_position()], \
				'right':[match.get_right_chromosome(), match.get_right_break_position()], \
				'id': match.uid, \
				'dataset':match.dataset_name }
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
	
	def set(self, \
			arg_left_chr, \
			arg_right_chr, \
			
			arg_left_pos, \
			arg_right_pos, \
			
			arg_left_strand, \
			arg_right_strand, \
			
			arg_uid, \
			arg_auto_set_acceptor_donor_direction
		):
		
		self.left_chr_str = arg_left_chr
		self.right_chr_str = arg_right_chr
		
		self.left_break_position = int(str(arg_left_pos).replace(",",""))
		self.right_break_position = int(str(arg_right_pos).replace(",",""))
		
		self.left_strand = arg_left_strand
		self.right_strand = arg_right_strand
		
		self.uid = arg_uid
		
		if (self.get_left_chromosome(False) > self.get_right_chromosome(False)) or \
		   ( \
				(self.get_left_chromosome(False) == self.get_right_chromosome(False)) and \
				(self.get_left_break_position() > self.get_right_break_position())
			):
			if self.acceptor_donor_direction == None and arg_auto_set_acceptor_donor_direction:
				self.acceptor_donor_direction = AD_DIRECTION_REVERSE
			
			self.swap_positions()
		else:
			if self.acceptor_donor_direction == None and\
			   arg_auto_set_acceptor_donor_direction:
				self.acceptor_donor_direction = AD_DIRECTION_FORWARD
	
	def find_strand_type(self,strand_type):
		if(strand_type in [STRAND_FORWARD, STRAND_REVERSE]):
			return strand_type
		
		if isinstance(strand_type, str):
			strand_type = strand_type.lower()
			if strand_type in ["f","+","forward","forwards","positive","5' -> 3'"]:
				return STRAND_FORWARD
			elif strand_type in ["b","-","r","backward","backwards","reverse","negative","3' -> 5'"]:
				return STRAND_REVERSE
			else:
				raise Exception("Unknown fusion strand: '"+strand_type+"'")
		
		return None
	
	def acceptor_donor_direction(self):
		return self.acceptor_donor_direction
	
	def cleanup_chr_name(self,chr_name):
		"""Given the large number of fusion genes, we remove all 'chr'
		prefixes because they add 6 bytes per fusion gene. They can be
		added again using the getter functions.
		"""
		
		chr_name = chr_name.strip()
		return chr_name[3:] if chr_name[0:3] == "chr" else chr_name
	
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
	
	def swap_positions(self):
		self.set( \
			self.right_chr_str, \
			self.left_chr_str, \
			self.get_right_break_position(), \
			self.get_left_break_position(), \
			self.right_strand, \
			self.left_strand, \
			self.uid, \
			self.acceptor_donor_direction != None
		)
	
	def annotate_genes_left(self,gene_names):
		self.annotated_genes_left = gene_names
		
	def annotate_genes_right(self,gene_names):
		self.annotated_genes_right = gene_names
	
	def get_annotated_genes_left(self,name_indexed):
		if(not name_indexed):
			if(not self.annotated_genes_left):
				return []
			else:
				return self.annotated_genes_left
		else:
			index = {}
			
			for gene in self.get_annotated_genes_left(False):
				gene_name = str(gene)
				if gene_name not in index:
					index[gene_name] = []
				index[gene_name].append(gene)
			
			return index
	
	def get_annotated_genes_right(self,name_indexed):
		if(not name_indexed):
			if(not self.annotated_genes_right):
				return []
			else:
				return self.annotated_genes_right
		else:
			index = {}
			for gene in self.get_annotated_genes_right(False):
				gene_name = str(gene)
				if gene_name not in index:
					index[gene_name] = []
				
				index[gene_name].append(gene)
			
			return index
	
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
	
	def get_left_strand(self):
		return self.left_strand
	
	def get_right_strand(self):
		return self.right_strand
	
	def has_annotated_genes(self):
		return self.annotated_genes_left and self.annotated_genes_right
	
	def show_me(self):
		print(self.__str__())
	
	def __str__(self):
		pos_left = self.get_left_position(True)
		pos_right = self.get_right_position(True)
		
		if(self.acceptor_donor_direction == AD_DIRECTION_FORWARD):
			acceptor_donor_direction = "->"
		elif(self.acceptor_donor_direction == AD_DIRECTION_REVERSE):
			acceptor_donor_direction = "<-"
		else:
			acceptor_donor_direction = "-"
		
		if(self.left_strand == STRAND_FORWARD):
			left_strand = "+"
		elif(self.left_strand == STRAND_REVERSE):
			left_strand = "-"
		else:
			left_strand = "?"
		
		if(self.right_strand == STRAND_FORWARD):
			right_strand = "+"
		elif(self.right_strand == STRAND_REVERSE):
			right_strand = "-"
		else:
			right_strand = "?"
		
		out = "Fusion '"+str(self.uid)+"' (from dataset '"+self.dataset_name+"'): " + self.get_left_chromosome(True)+":"+str(pos_left[1])+"("+left_strand+")" + acceptor_donor_direction + self.get_right_chromosome(True)+":"+str(pos_right[1])+"("+right_strand+")"
		if(self.get_annotated_genes_left(False)):
			out += "\n - annotated genes left:  "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_left(False)])
		if(self.get_annotated_genes_right(False)):
			out += "\n - annotated genes right: "+", ".join([str(gene_name) for gene_name in self.get_annotated_genes_right(False)])
		
		return out
