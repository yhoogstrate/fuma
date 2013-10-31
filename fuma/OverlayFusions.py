#!/usr/bin/env python



from ReadCGhighConfidenceJunctionsBeta import ReadCGhighConfidenceJunctionsBeta
from ReadTophatFusion import ReadTophatFusion
from ReadDefuse import ReadDefuse
from ParseBED import ParseBED
from CompareFusionsGTFOverlay import CompareFusionsGTFOverlay



import os.path



class OverlayFusions:
	def __init__(self):
		self.datasets = []
		self.matches = False
	
	def add_dataset(self,dataset):
		self.datasets.append(dataset)
	
	def create_matrix(self,n):
		matrix = []
		for dataset in self.datasets:
			matrix.append([])
		
		n = len(self.datasets)-1
		while(n >= 1):
			for i in range(n):
				matrix[n].append(i)
			n -= 1
		
		return matrix
	
	def overlay_fusions(self):
		matrix = self.create_matrix(len(self.datasets))
		
		for i in range(len(matrix)):
			for j in range(len(matrix[i])):
				matches = self.compare_datasets(self.datasets[i],self.datasets[j])
				matrix[i][j] = matches
		
		self.matches = matrix
		return self.matches
	
	def compare_datasets(self,dataset_1,dataset_2):
		types = sorted([dataset_1.get_type(),dataset_2.get_type()])
		
		if(types[0] == "DNA" and types[1] == "DNA"):
			#comparison = CompareFusionsPositionOverlay(dataset_1,dataset_2,insertsize)
			pass
		elif(types[0] == "DNA" and types[1] == "RNA"):
			comparison = CompareFusionsGTFOverlay(dataset_1,dataset_2,self.gene_annotation)
		elif(types[0] == "RNA" and types[1] == "RNA"):
			comparison = CompareFusionsGTFOverlay(dataset_1,dataset_2,self.gene_annotation)
		
		matches = comparison.compare()
		return matches
	
	def set_annotation(self,arg_gene_annotation):
		self.gene_annotation = arg_gene_annotation
	
	def export(self,filename_prefix="",join="_vs._",suffix=".txt"):
		for i in range(len(self.matches)):
			for j in range(len(self.matches[i])):
				filename = os.path.basename(filename_prefix+self.datasets[i].name)+join+os.path.basename(self.datasets[j].name)+suffix
				print "exporting: "+os.path.basename(filename)
				fh = open(filename,"w")
				
				fh.write("["+self.datasets[i].name+"]-position\t")
				fh.write("["+self.datasets[i].name+"]-left-junction-associated-genes\t")
				fh.write("["+self.datasets[i].name+"]-right-junction-associated-genes\t")
				fh.write("["+self.datasets[j].name+"]-position\t")
				fh.write("["+self.datasets[j].name+"]-left-junction-associated-genes\t")
				fh.write("["+self.datasets[j].name+"]-right-junction-associated-genes\n")
				
				for match in self.matches[i][j]:
					fh.write(str(match[0].get_left_chromosome())+":"+str(match[0].get_left_break_position())+"-"+str(match[0].get_right_chromosome())+":"+str(match[0].get_right_break_position())+"\t")
					fh.write(";".join(match[0].get_annotated_genes_left())+"\t")
					fh.write(";".join(match[0].get_annotated_genes_right())+"\t")
					fh.write(str(match[1].get_left_chromosome())+":"+str(match[1].get_left_break_position())+"-"+str(match[1].get_right_chromosome())+":"+str(match[1].get_right_break_position())+"\t")
					fh.write(";".join(match[1].get_annotated_genes_left())+"\t")
					fh.write(";".join(match[1].get_annotated_genes_right())+"\n")
				fh.close()
