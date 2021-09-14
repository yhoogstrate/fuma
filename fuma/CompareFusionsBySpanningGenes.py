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

import logging

from .Fusion import Fusion
from .FusionDetectionExperiment import FusionDetectionExperiment


class CompareFusionsBySpanningGenes:
	logger = logging.getLogger("FuMa::CompareFusionsBySpanningGenes")
	
	def __init__(self,experiment_1,experiment_2,args):
		self.experiment_1 = experiment_1
		self.experiment_2 = experiment_2
		
		self.args = args
	
	def find_overlap(self):
		self.logger.debug("Comparing: '"+self.experiment_1.name+"' with '"+self.experiment_2.name + "'" + " - using '"+self.args.matching_method+"'-based matching")
		overlap_between_experiments = FusionDetectionExperiment(self.experiment_1.name+"_vs._"+self.experiment_2.name)
		
		if(self.experiment_1.genes_spanning_left_junction and self.experiment_2.genes_spanning_left_junction and self.experiment_1.genes_spanning_right_junction and self.experiment_2.genes_spanning_right_junction):
			overlap_between_experiments.genes_spanning_left_junction = list(set(self.experiment_1.genes_spanning_left_junction+self.experiment_2.genes_spanning_left_junction))
			overlap_between_experiments.genes_spanning_right_junction = list(set(self.experiment_1.genes_spanning_right_junction+self.experiment_2.genes_spanning_right_junction))
			
			matches_exp_1 = set()
			matches_exp_2 = set()
			
			for chromosome_left in self.experiment_1.index.items():
				for chromosome_right in chromosome_left[1].items():
					for fusion_1 in chromosome_right[1]:
						
						if chromosome_left[0] in self.experiment_2.index and chromosome_right[0] in self.experiment_2.index[chromosome_left[0]]:
							for fusion_2 in self.experiment_2.index[chromosome_left[0]][chromosome_right[0]]:
								
								## Do the gene-name comparison
								#if(self.args.matching_method == 'egm'):
								#	match = self.match_fusions_egm(fusion_1,fusion_2,False)
								#else:
								match = self.match_fusions(fusion_1,fusion_2,False)
								
								if(match):
									match.matches = fusion_1.matches | fusion_2.matches
									
									matches_exp_1.add(fusion_1)
									matches_exp_2.add(fusion_2)
									
									fusion_1.matched_datasets[fusion_2.dataset_name] = True
									fusion_2.matched_datasets[fusion_1.dataset_name] = True
									
									overlap_between_experiments.add_fusion(match)
			
			overlap_between_experiments.remove_duplicates(self.args)
			
			return [overlap_between_experiments,len(matches_exp_1),len(matches_exp_2),(matches_exp_1 | matches_exp_2)]
		else:
			self.logger.warning("No gene annotation reference found")
	
	
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
	
	
	def match_fusions(self,fusion_1,fusion_2,allow_empty = True):
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
			
			if((fusion_1.annotated_genes_left and fusion_1.annotated_genes_right and fusion_2.annotated_genes_left and fusion_2.annotated_genes_right)):
				# Check if all of the smallest are in the largest;
				# if you do it otherwise you don't know if all from the smallest are also in the largest
				
				fusion_1_annotated_genes_left =  fusion_1.get_annotated_genes_left(True)
				fusion_1_annotated_genes_right = fusion_1.get_annotated_genes_right(True)
				
				fusion_2_annotated_genes_left =  fusion_2.get_annotated_genes_left(True)
				fusion_2_annotated_genes_right = fusion_2.get_annotated_genes_right(True)
				
				if(self.args.matching_method == 'overlap'):
					matches_left  = self.match_overlap( set(fusion_1_annotated_genes_left.keys()), set(fusion_2_annotated_genes_left.keys()) )
					matches_right = self.match_overlap( set(fusion_1_annotated_genes_right.keys()), set(fusion_2_annotated_genes_right.keys()) )
				elif(self.args.matching_method == 'egm'):
					matches_left  = self.match_egm( set(fusion_1_annotated_genes_left.keys()), set(fusion_2_annotated_genes_left.keys()) )
					matches_right = self.match_egm( set(fusion_1_annotated_genes_right.keys()), set(fusion_2_annotated_genes_right.keys()) )
				else:
					matches_left  = self.match_sets( set(fusion_1_annotated_genes_left.keys()), set(fusion_2_annotated_genes_left.keys()) )
					matches_right = self.match_sets( set(fusion_1_annotated_genes_right.keys()), set(fusion_2_annotated_genes_right.keys()) )
				
				# Do we allow empty matches as empty results or 2x empty input? >> if the latter, the if should be in the beginning of the function
				if(matches_left and matches_right and \
						(allow_empty or \
						(len(fusion_1.annotated_genes_left) > 0 and \
						len(fusion_1.annotated_genes_right) > 0 and \
						len(fusion_2.annotated_genes_left) > 0 and \
						len(fusion_2.annotated_genes_right) > 0)) \
					):
					
					fusion_merged = Fusion( \
						fusion_1.get_left_chromosome(), \
						fusion_1.get_right_chromosome(), \
						fusion_1.get_left_break_position(), \
						fusion_1.get_right_break_position(), \
						fusion_1.left_strand, \
						fusion_1.right_strand, \
						fusion_1.dataset_name+"_vs._"+fusion_2.dataset_name, \
						"", \
						(fusion_1.acceptor_donor_direction != None and fusion_2.acceptor_donor_direction != None)
					)
					
					# Fancy oneliners: create a list of all Gene objects based on the gene names in matches_left
					#print [item for sublist in matches_left for item in fusion_1_annotated_genes_left[sublist]+fusion_2_annotated_genes_left[sublist]]
					
					fusion_merged.annotate_genes_left(list(set([item for sublist in matches_left for item in ((fusion_1_annotated_genes_left[sublist]) if sublist in fusion_1_annotated_genes_left else [])+((fusion_2_annotated_genes_left[sublist]) if sublist in fusion_2_annotated_genes_left else [])])))
					fusion_merged.annotate_genes_right(list(set([item for sublist in matches_right for item in ((fusion_1_annotated_genes_right[sublist]) if sublist in fusion_1_annotated_genes_right else [])+((fusion_2_annotated_genes_right[sublist]) if sublist in fusion_2_annotated_genes_right else [])])))
					
					#@todo Check whether keeping the references to the original fusion objects is much more intensive or not - if not, use it instead
					#   otherwise, make a Location() object and use it and save the reference in the fusion class
					#for location in fusion_1.locations()+fusion_2.locations():
					#	fusion_merged.add_location(location)
					
					# If one fusion is (A,B) and the other (B,A), the directions are opposite
					# Therefore not the direction of fusion_1 should be chosen, but it should be set to "None" / unknown
					acceptor_donor_directions = set([fusion_1.acceptor_donor_direction,fusion_2.acceptor_donor_direction])
					if(len(acceptor_donor_directions) != 1):
						fusion_merged.acceptor_donor_direction = None
					else:
						fusion_merged.acceptor_donor_direction = list(acceptor_donor_directions)[0]
					
					# If one fusion's left strand is (+) and the other is (-) the strand should be unknown, similarly for the right
					left_strands = set([fusion_1.left_strand,fusion_2.left_strand])
					if(len(left_strands) != 1):
						fusion_merged.left_strand = None
					else:
						fusion_merged.left_strand = list(left_strands)[0]
						
					right_strands = set([fusion_1.right_strand,fusion_2.right_strand])
					if(len(right_strands) != 1):
						fusion_merged.right_strand = None
					else:
						fusion_merged.right_strand = list(right_strands)[0]
					
					return fusion_merged
				else:
					return False
			else:
				return False
		else:
			return False
	
	def match_sets(self,superset,subset):								#https://docs.python.org/2/library/sets.html
		if(len(subset) > len(superset)):
			return self.match_sets(subset,superset)						# Gene names have to be provided as sets
		elif(subset.issubset(superset)):
			return subset
		else:
			return None
	
	def match_egm(self,set1,set2):
		if set1 == set2:
			return set1
		else:
			return None
	
	def match_overlap(self,set1,set2):									#https://docs.python.org/2/library/sets.html
		overlap = set1.intersection(set2)
		if(not overlap):
			return None
		else:
			return overlap
	
	""" 
	#This type of matching increases the sets after multiple iterations
	def match_sets_return_superset(self,superset,subset):
		if(subset == superset):
			return superset												#(subset | superset)
		else:
			return False
	"""
