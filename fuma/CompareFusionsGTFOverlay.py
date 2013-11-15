#!/usr/bin/env python



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class CompareFusionsGTFOverlay:
	def __init__(self,dataset_1,dataset_2):
		self.dataset_1 = dataset_1
		self.dataset_2 = dataset_2
	
	def compare(self):
		print " - Comparing: "+self.dataset_1.name+" - "+self.dataset_2.name
		
		dataset = HighThroughputFusionDetectionExperiment(self.dataset_1.name+"_vs._"+self.dataset_2.name,"RNA")
		
		matches = []
		
		for chromosome_1 in self.dataset_1.get_fusions():
			for fusion_1 in chromosome_1["fusions"]:
				left_genes = []
				right_genes = []
				
				found = False
				
				for chromosome_2 in self.dataset_2.get_fusions():
					if(chromosome_1["name"] == chromosome_2["name"]):
						for fusion_2 in chromosome_2["fusions"]:
							match = self.match_fusions(fusion_1,fusion_2,False)
							if(match):
								matches.append([fusion_1,fusion_2])
								dataset.add_fusion(match)
		
		dataset = self.remove_duplicates(dataset)
		
		return [matches,dataset]
	
	def remove_duplicates(self,arg_dataset):
		"""
		- First create a table of those that overlap
		- Then create merged entries based on the overlap matrix
		"""
		
		dataset = HighThroughputFusionDetectionExperiment(arg_dataset.name,arg_dataset.type)
		
		for chromosome in arg_dataset.get_fusions():
			fusions = chromosome["fusions"]
			for i in range(len(fusions)):
				fusion_1 = fusions[i]
				if(fusion_1):
					for j in range(i+1,len(fusions)):
						fusion_2 = fusions[j]
						if(fusion_2):
							match = self.match_fusions(fusion_1,fusion_2,True)
							if(match):
								fusion_1 = match
								fusions[j] = False
					dataset.add_fusion(fusion_1)
		return dataset
	
	def match_fusions(self,fusion_1,fusion_2,allow_empty = True):
		"""Matches whether two fusion objects are the same prediction
				# fusion_1 <=> fusion_2; for both left and right position:
				# [a,b,c] == [a,b,c]     ->    [a,b,c]
				# [a,b,c] == [a,b]       ->    [a,b,c]
				# 
				
				# [a,b,c,d] != [a,b,e]
				#	BECAUSE: [a,b] can not be located in C, never
				#
		"""
		
		matches_left = False
		matches_right = False
		
		if(type(fusion_1.get_annotated_genes_left()) == type([]) and type(fusion_1.get_annotated_genes_right()) == type([]) and type(fusion_2.get_annotated_genes_left()) == type([]) and type(fusion_2.get_annotated_genes_right()) == type([])):
			matches_left = True
			matches_right = True
			
			# Check if all of the smallest are in the largest; if you do it otherwise you don't know if all from the smallest are also in the largest
			if(len(fusion_1.get_annotated_genes_left()) <  len(fusion_2.get_annotated_genes_left())):
				fusions_left_small = fusion_1.get_annotated_genes_left()
				fusions_left_large = fusion_2.get_annotated_genes_left()
			else:
				fusions_left_small = fusion_2.get_annotated_genes_left()
				fusions_left_large = fusion_1.get_annotated_genes_left()
			
			if(len(fusion_1.get_annotated_genes_right()) <  len(fusion_2.get_annotated_genes_right())):
				fusions_right_small = fusion_1.get_annotated_genes_right()
				fusions_right_large = fusion_2.get_annotated_genes_right()
			else:
				fusions_right_small = fusion_2.get_annotated_genes_right()
				fusions_right_large = fusion_1.get_annotated_genes_right()
			
			for gene_1 in fusions_left_small:
				if(gene_1 not in fusions_left_large):
					matches_left = False
			
			for gene_1 in fusions_right_small:
				if(gene_1 not in fusions_right_large):
					matches_right = False
		
		if(matches_left and matches_right):
			if(allow_empty or (len(fusions_left_large) > 0 and len(fusions_right_large) > 0)):
				print fusion_1.get_left_chromosome(),fusion_1.get_left_chromosome(True)
				
				fusion_merged = Fusion(fusion_1.get_left_chromosome(),fusion_1.get_right_chromosome(),fusion_1.get_left_break_position(),fusion_1.get_right_break_position(),fusion_1.get_sequence(),fusion_1.get_transition_sequence(),fusion_1.get_left_strand(),fusion_1.get_right_strand(),fusion_1.get_dataset_name())
				fusion_merged.annotate_genes_left(fusions_left_large)
				fusion_merged.annotate_genes_right(fusions_right_large)
				
				for key in fusion_1.locations.keys():
					fusion_merged.set_locations(key,fusion_1.locations[key])
				for key in fusion_2.locations.keys():
					fusion_merged.set_locations(key,fusion_2.locations[key])
				
				fusion_merged.dataset_name = fusion_1.get_dataset_name()+"_vs._"+fusion_2.get_dataset_name()
				
				return fusion_merged
			else:
				return False
		else:
			return False
