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

from Readers import *

from ParseBED import ParseBED
from FusionDetectionExperiment import FusionDetectionExperiment
from Triangle import Triangle
from MergedFusion import MergedFusion


import os.path,sys,itertools


class ComparisonMatrix:
	logger = logging.getLogger("FuMa::ComparisonMatrix")
	
	def __init__(self,args):
		self.experiments = []
		self.args = args
	
	def add_experiment(self,arg_experiment):
		if not isinstance(arg_experiment, FusionDetectionExperiment):
			raise Exception("MergedFusion objects can only be expanded with Fusion objects")
		else:
			if arg_experiment in self.experiments:
				raise Exception("MergedFusion is updated with one that it already contains")
			else:
				self.experiments.append(arg_experiment)
	
	def overlay_fusions(self):
		fusions = []
		
		for experiment in self.experiments:
			for fusion in experiment:
				fusions.append(fusion)
		
		n = len(fusions)
		
		"""
		iter1: 
		
		0,7 | 1,7 2,7 | 3,7 4,7 5,7 | 6,7 7,7
		
		iter2:
		0,6 | 1,6 2,6 | 3,6 4,6 5,6 | 6,6
		
		"""
		for y in range(len(fusions)-1,0,-1):
			for x in range(y+1):
				# If they do not belong to the same dataset - i.e. no duplication removal
				# And if they are the same MergedFusion gene
				if self.map_i_to_exp_id(x) != self.map_i_to_exp_id(y) and fusions[x] != fusions[y]:
					## do actual comparison
					#comparison, additional_replacements = self.match_fusions(fusions[x], fusions[y])
					comparison = self.match_fusions(fusions[x], fusions[y])
					
					if comparison != False:
						fusions[x] = comparison
						fusions[y] = comparison
						
						# In case a MergedFusion is identical to another MergedFusion, one of the objects needs te be replaced with the other
						#if additional_replacements:
						#	for z in self.get_merged_fusion_occurances(fusions, additional_replacements):
						#		fusions[z] = comparison
						#	del(additional_replacements)
		
		self.export_list(fusions)
	
	def export_list(self,fusions):
		todo = len(fusions)
		k = 1
		while todo > 0:
			print "k=",k
			for i in range(len(fusions)):
				fusion = fusions[i]
				if fusion != None:
					if k == 1 and isinstance(fusion, Fusion):
						print fusion
						fusions[i] = None
						todo -= 1
					else:
						if len(fusions[i]) == k:
							print fusion
							for j in range(len(fusions)):
								if fusions[j] == fusion:
									fusions[j] = None
									todo -= 1
							del(fusion)
			k += 1
		print "done"
	
	def map_i_to_exp_id(self,i):
		"""
		Maps the i-th position back to the id of the corresponding experiment
		"""
		cumulative = 0
		j = -1
		while i >= cumulative:
			j += 1
			cumulative += len(self.experiments[j])
		return j
	
	def get_merged_fusion_occurances(self, fusions, additional_replacements):
		for i in range(len(fusions)):
			if fusions[i] == additional_replacements:
				yield i
	
	def __len__(self):
		return len(self.experiments)
	
	def match_fusions(self,fusion_1,fusion_2):
		"""Matches whether two fusion objects are the same prediction
				# fusion_1 <=> fusion_2; for both left and right position:
				# [a,b,c] == [a,b,c]     ->    [a,b,c]
				# [a,b,c] == [a,b]       ->    [a,b,c]
				# 
				
				# [a,b,c,d] != [a,b,e]
				#	BECAUSE: [a,b] can not be located in C, never
				#
				# (not is_empty(a)) and subset(a,b) or subset(b,a)
		"""
		
		# First check whether the strands match, if strand-specific-matching is enabled:
		if(self.match_fusion_gene_strands(fusion_1,fusion_2) and self.match_acceptor_donor_direction(fusion_1,fusion_2)):
			if fusion_1.has_annotated_genes() and fusion_2.has_annotated_genes():
				# Check if all of the smallest are in the largest;
				# if you do it otherwise you don't know if all from the smallest are also in the largest
				
				if(self.args.matching_method == 'overlap'):
					matches_left  = self.match_overlap(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
					matches_right = self.match_overlap(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
				elif(self.args.matching_method == 'egm'):
					matches_left  = self.match_egm(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
					matches_right = self.match_egm(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
				else:
					matches_left  = self.match_sets(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
					matches_right = self.match_sets(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
				
				# Do we allow empty matches as empty results or 2x empty input? >> if the latter, the if should be in the beginning of the function
				if matches_left and matches_right:
					#replace_merged_fusions = None
					if isinstance(fusion_1, MergedFusion) and isinstance(fusion_2, MergedFusion):
						raise Exception("If (A & B) == (C & D), (A & B & C) should have matched before..")
						#merged_fusion = fusion_1
						#merged_fusion.merge(fusion_2)
						#replace_merged_fusions = fusion_2
					elif isinstance(fusion_1, MergedFusion) and isinstance(fusion_2, Fusion):
						merged_fusion = fusion_1
						merged_fusion.add_fusion(fusion_2)
					elif isinstance(fusion_1, Fusion) and isinstance(fusion_2, MergedFusion):
						merged_fusion = fusion_2
						merged_fusion.add_fusion(fusion_1)
					elif isinstance(fusion_1, Fusion) and isinstance(fusion_2, Fusion):
						merged_fusion = MergedFusion()
						merged_fusion.add_fusion(fusion_1)
						merged_fusion.add_fusion(fusion_2)
					else:
						raise Exception("Something went wrong with the object types")
					
					# This has to be pre-cached and can not be determined on the fly by a functions,
					# because it requires the type of matching. If you would allow for functions, you could 
					# end up with overlap and egm and subset based matching mixed up.
					merged_fusion.annotated_genes_left = matches_left
					merged_fusion.annotated_genes_right = matches_right
					
					return merged_fusion#, replace_merged_fusions
				else:
					return False#, None
			else:
				return False#, None
		else:
			return False#, None
	
	def match_fusion_gene_strands(self,fusion_1,fusion_2):
		if not self.args.strand_specific_matching:
			return True
		else:
			if fusion_1.left_strand == None or fusion_1.right_strand == None or fusion_2.left_strand == None or fusion_2.right_strand == None:
				raise Exception("A fusion gene without an annotated strand was used for strand-specific-matching.\n\n"+fusion_1.__str__()+"\n"+fusion_2.__str__())
			else:
				return fusion_1.left_strand == fusion_2.left_strand and fusion_1.right_strand == fusion_2.right_strand
	
	def match_acceptor_donor_direction(self,fusion_1,fusion_2):
		if(not self.args.acceptor_donor_order_specific_matching):
			return True
		elif(fusion_1.acceptor_donor_direction == None or fusion_2.acceptor_donor_direction == None):
			raise Exception("A fusion gene without an annotated acceptor-donor direction was used for acceptor-donor-order-specific-matching.\n\n"+fusion_1.__str__()+"\n"+fusion_2.__str__())
		else:
			return (fusion_1.acceptor_donor_direction == fusion_2.acceptor_donor_direction)
	
	def match_sets(self,superset,subset):								#https://docs.python.org/2/library/sets.html
		subset_s = set([str(gene) for gene in subset])
		superset_s = set([str(gene) for gene in superset])
		
		if(len(subset_s) > len(superset_s)):
			return self.match_sets(subset,superset)						# Gene names have to be provided as sets
		elif(subset_s.issubset(superset_s)):
			return subset
		else:
			return None
