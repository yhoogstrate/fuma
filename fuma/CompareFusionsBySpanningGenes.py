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

from Fusion import Fusion
from FusionDetectionExperiment import FusionDetectionExperiment


class CompareFusionsBySpanningGenes:
	def __init__(self,experiment_1,experiment_2):
		self.experiment_1 = experiment_1
		self.experiment_2 = experiment_2
	
	def find_overlap(self):
		print " - Comparing: "+self.experiment_1.name+" - "+self.experiment_2.name
		
		overlap_between_experiments = FusionDetectionExperiment(self.experiment_1.name+"_vs._"+self.experiment_2.name,"RNA")
		
		if(self.experiment_1.genes_spanning_left_junction and self.experiment_2.genes_spanning_left_junction and self.experiment_1.genes_spanning_right_junction and self.experiment_2.genes_spanning_right_junction):
			overlap_between_experiments.genes_spanning_left_junction = list(set(self.experiment_1.genes_spanning_left_junction+self.experiment_2.genes_spanning_left_junction))
			overlap_between_experiments.genes_spanning_right_junction = list(set(self.experiment_1.genes_spanning_right_junction+self.experiment_2.genes_spanning_right_junction))
			
			matches_exp_1 = []
			matches_exp_2 = []
			
			for chromosome_left in self.experiment_1.index.items():
				for chromosome_right in chromosome_left[1].items():
					for fusion_1 in chromosome_right[1]:
						
						if(self.experiment_2.index.has_key(chromosome_left[0]) and self.experiment_2.index[chromosome_left[0]].has_key(chromosome_right[0])):
							for fusion_2 in self.experiment_2.index[chromosome_left[0]][chromosome_right[0]]:
								match = self.match_fusions(fusion_1,fusion_2,False)
								if(match):
									matches_exp_1.append(fusion_1)
									matches_exp_2.append(fusion_2)
									overlap_between_experiments.add_fusion(match)
									found = True
			
			overlap_between_experiments.remove_duplicates("by-gene-names")
			
			return [overlap_between_experiments,len(list(set(matches_exp_1))),len(list(set(matches_exp_2)))]
		else:
			print "No gene annotation reference found"
	
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
		
		if((fusion_1.annotated_genes_left and fusion_1.annotated_genes_right and fusion_2.annotated_genes_left and fusion_2.annotated_genes_right)):
			
			# Check if all of the smallest are in the largest;
			# if you do it otherwise you don't know if all from the smallest are also in the largest
			
			fusion_1_annotated_genes_left =  fusion_1.get_annotated_genes_left(True)
			fusion_1_annotated_genes_right = fusion_1.get_annotated_genes_right(True)
			
			fusion_2_annotated_genes_left =  fusion_2.get_annotated_genes_left(True)
			fusion_2_annotated_genes_right = fusion_2.get_annotated_genes_right(True)
			
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
				
				fusion_merged = Fusion(fusion_1.get_left_chromosome(),fusion_1.get_right_chromosome(),fusion_1.get_left_break_position(),fusion_1.get_right_break_position(),fusion_1.get_sequence(),fusion_1.get_transition_sequence(),fusion_1.get_left_strand(),fusion_1.get_right_strand(),fusion_1.get_dataset_name()+"_vs._"+fusion_2.get_dataset_name())
				
				# Awesome oneliner: create a list of all Gene objects based on the gene names in matches_left
				#print [item for sublist in matches_left for item in fusion_1_annotated_genes_left[sublist]+fusion_2_annotated_genes_left[sublist]]
				
				fusion_merged.annotate_genes_left(list(set([item for sublist in matches_left for item in ((fusion_1_annotated_genes_left[sublist]) if fusion_1_annotated_genes_left.has_key(sublist) else [])+((fusion_2_annotated_genes_left[sublist]) if fusion_2_annotated_genes_left.has_key(sublist) else [])])))
				
				fusion_merged.annotate_genes_right(list(set([item for sublist in matches_right for item in ((fusion_1_annotated_genes_right[sublist]) if fusion_1_annotated_genes_right.has_key(sublist) else [])+((fusion_2_annotated_genes_right[sublist]) if fusion_2_annotated_genes_right.has_key(sublist) else [])])))
				
				#@todo Check whether keeping the references to the original fusion objects is much more intensive or not - if not, use it instead
				#   otherwise, make a Location() object and use it and save the reference in the fusion class
				for location in fusion_1.locations+fusion_2.locations:
					fusion_merged.add_location(location)
				
				return fusion_merged
			else:
				return False
		else:
			return False
	
	def match_sets(self,superset,subset):								#https://docs.python.org/2/library/sets.html
		#return self.match_equal_sets(superset,subset)
		
		if(len(subset) > len(superset)):
			return self.match_sets(subset,superset)						# Gene names have to be provided as sets?!
		else:
			if(subset.issubset(superset)):
				return subset
			else:
				return False
	
	def match_equal_sets(self,superset,subset):
		if(subset == superset):
			return subset
		else:
			return False
	
	def match_sets_return_superset(self,superset,subset):
		if(subset == superset):
			return superset												#(subset | superset)
		else:
			return False
