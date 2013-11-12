#!/usr/bin/env python



class CompareFusionsGTFOverlay:
	def __init__(self,dataset_1,dataset_2):
		self.dataset_1 = dataset_1
		self.dataset_2 = dataset_2
	
	def compare(self):
		print " - Comparing: "+self.dataset_1.name+" - "+self.dataset_2.name
		
		matches = []
		#aggregated_match_object = Highthroughput...()
		
		for chromosome_1 in self.dataset_1.get_fusions():
			for fusion_1 in chromosome_1["fusions"]:
				#for fusion_1 in position_1:
				left_genes = []
				right_genes = []
				
				found = False
				
				for chromosome_2 in self.dataset_2.get_fusions():
					if(chromosome_1["name"] == chromosome_2["name"]):
						for fusion_2 in chromosome_2["fusions"]:
							match = self.match_fusions(fusion_1,fusion_2)
							if(match):
								matches.append([fusion_1,fusion_2])
								
								#fusion = Fusion()
								#fusion.annotate_left_genes(match[["left"])
								#fusion.annotate_right_genes(match[["right"])
								
								#aggregated_match_object.add_fusion(fusion)
				
				# fusion_1 <=> fusion_2
				# [a,b,c] == [a,b,c]     ->    [a,b,c]
				# [a,b,c] == [a,b]     ->      ?   [a,b] ? [a,b,c] ? [a,b]&[a,b,c]
				
				# [a,b,c] != [a,b]
				#	BECAUSE: [a,b] can not be located in C, never
				#
		
		# Matches as Table
		# Matches as "New fusion dataset" <- essential for comparing multiple with each other
		
		return matches
		#return [matches,aggregated_match_object]
	
	def match_fusions(self,fusion_1,fusion_2):
		matches_left = []
		matches_right = []
		
		if(type(fusion_1.get_annotated_genes_left()) == type([]) and type(fusion_1.get_annotated_genes_right()) == type([]) and type(fusion_2.get_annotated_genes_left()) == type([]) and type(fusion_2.get_annotated_genes_right()) == type([])):
			for gene_1 in fusion_1.get_annotated_genes_left():
				for gene_2 in fusion_2.get_annotated_genes_left():
					if(gene_1 == gene_2):
						matches_left.append(gene_1)
			
			for gene_1 in fusion_1.get_annotated_genes_right():
				for gene_2 in fusion_2.get_annotated_genes_right():
					if(gene_1 == gene_2):
						matches_right.append(gene_1)
		
		if(len(matches_left) > 0 and len(matches_right) > 0):
			return {"left":matches_left,"right":matches_right}
		else:
			return False
