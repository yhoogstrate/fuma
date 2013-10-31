#!/usr/bin/env python



class CompareFusionsGTFOverlay:
	def __init__(self,dataset_1,dataset_2,arg_gene_annotation):
		self.dataset_1 = dataset_1
		self.dataset_2 = dataset_2
		self.gene_annotation = arg_gene_annotation
	
	def compare(self):
		"""
		@ todo  INDEX genes by name-sorting?
		"""
		print " - Comparing: "+self.dataset_1.name+" - "+self.dataset_2.name
		
		self.gene_annotation.convert_left_locations_to_genes(self.dataset_1)
		self.gene_annotation.convert_left_locations_to_genes(self.dataset_2)
		
		self.gene_annotation.convert_right_locations_to_genes(self.dataset_1)
		self.gene_annotation.convert_right_locations_to_genes(self.dataset_2)
		
		
		matches = []
		
		for chromosome_1 in self.dataset_1.get_fusions():
			for fusion_1 in chromosome_1["fusions"]:
				#for fusion_1 in position_1:
				found = False
				
				for chromosome_2 in self.dataset_2.get_fusions():
					if(chromosome_1["name"] == chromosome_2["name"]):
						for fusion_2 in chromosome_2["fusions"]:
							#for fusion_2 in position_2:
							"""
								if(fusion_1.get_left_chromosome(False) == 17 and fusion_1.get_right_chromosome(False) == 20):
									if(fusion_2.get_left_chromosome(False) == 17 and fusion_2.get_right_chromosome(False) == 20):
										print "1"
										fusion_1.show_me()
										print "2"
										fusion_2.show_me()
										print "match: ",self.match_fusions(fusion_1,fusion_2)
										print
										print
							"""
							if(self.match_fusions(fusion_1,fusion_2)):
								matches.append([fusion_1,fusion_2])
		
		return matches
	
	def match_fusions(self,fusion_1,fusion_2):
		match_left = False
		match_right = False
		
		if(type(fusion_1.get_annotated_genes_left()) == type([]) and type(fusion_1.get_annotated_genes_right()) == type([]) and type(fusion_2.get_annotated_genes_left()) == type([]) and type(fusion_2.get_annotated_genes_right()) == type([])):
			for gene_1 in fusion_1.get_annotated_genes_left():
				for gene_2 in fusion_2.get_annotated_genes_left():
					if(gene_1 == gene_2):
						match_left = True
			
			for gene_1 in fusion_1.get_annotated_genes_right():
				for gene_2 in fusion_2.get_annotated_genes_right():
					if(gene_1 == gene_2):
						match_right = True
		
		return (match_left and match_right)
