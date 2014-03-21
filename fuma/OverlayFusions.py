#!/usr/bin/env python



from Readers import *

from ParseBED import ParseBED
#from CompareFusionsDistance import CompareFusionsDistance
from CompareFusionsGTFOverlay import CompareFusionsGTFOverlay



import os.path
import itertools



class OverlayFusions:
	def __init__(self):
		self.datasets = []
		self.matches = False
		self.distances = False
	
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
	
	def create_keys(self,combination):
		keys = [combination[:-1],combination[-1:],combination]
		
		for i in range(len(keys)):
			for j in range(len(keys[i])):
				keys[i][j] = str(keys[i][j])
			
			keys[i] = ".".join(keys[i])
		
		return keys
	
	def overlay_fusions(self):
		n = len(self.datasets)
		
		self.scores_tmp = {}
		self.matrix_tmp = {}
		
		for i in range(len(self.datasets)):
			self.matrix_tmp[str(i+1)] = self.datasets[i]
		
		comparisons = self.find_combination_table(n)
		
		for r in comparisons:
			for c in r:
				keys = self.create_keys(c)
				
				matches = self.compare_datasets(self.matrix_tmp[keys[0]],self.matrix_tmp[keys[1]])
				
				self.matrix_tmp[keys[2]] = matches[1]
				self.scores_tmp[keys[2]] = matches[0]
		
		return self.scores_tmp
	
	def find_combination_table(self,n):
		in_list = range(1,n+1)
		
		table = []
		
		for r in range(2,len(in_list)+1):
			table_i_tmp = itertools.combinations(in_list,r)
			table_i = []
			for i in table_i_tmp:
				table_i.append(list(i))
			table.append(table_i)
		
		return table
	
	"""
	def find_distances(self):
		self.distances = self.create_matrix(len(self.datasets))
		
		for i in range(len(self.distances)):
			for j in range(len(self.distances[i])):
				o = CompareFusionsDistance(self.datasets[i],self.datasets[j])
				self.distances[i][j] = o.compare()
	"""
	
	def compare_datasets(self,dataset_1,dataset_2):
		types = sorted([dataset_1.get_type(),dataset_2.get_type()])
		
		if(types[0] == "DNA" and types[1] == "DNA"):
			comparison = CompareFusionsGTFOverlay(dataset_1,dataset_2)
			#comparison = CompareFusionsPositionOverlay(dataset_1,dataset_2,insertsize)
			#pass
		elif(types[0] == "DNA" and types[1] == "RNA"):
			comparison = CompareFusionsGTFOverlay(dataset_1,dataset_2)
		elif(types[0] == "RNA" and types[1] == "RNA"):
			comparison = CompareFusionsGTFOverlay(dataset_1,dataset_2)
		
		matches = comparison.compare()
		return matches
	
	def set_annotation(self,arg_gene_annotation):
		self.gene_annotation = arg_gene_annotation
	
	def export(self,filename_prefix="",join="_vs._",suffix=".txt"):
		for i in range(len(self.matches)):
			for j in range(len(self.matches[i])):
				filename = filename_prefix+self.datasets[i].name+join+os.path.basename(self.datasets[j].name)+suffix
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
	
	def export2(self,filename_prefix="",suffix=".txt"):
		for key in self.scores_tmp.keys():
			name = self.matrix_tmp[key].name
			filename = filename_prefix+os.path.basename(name)+suffix
			print "exporting: "+os.path.basename(filename)
			fh = open(filename,"w")
			
			datasets = {}
			
			for chromosome in self.matrix_tmp[key].get_fusions():
				for fusion in chromosome["fusions"]:
					for location in fusion.locations:
						datasets[location["dataset"]] = True
			datasets = sorted(datasets.keys())
			
			fh.write("genes-left-junction\tgenes-right-junction\t"+"\t".join(datasets)+"\n")
			
			for chromosome in self.matrix_tmp[key].get_fusions():
				for fusion in chromosome["fusions"]:
					fh.write(";".join(fusion.get_annotated_genes_left())+"\t")
					fh.write(";".join(fusion.get_annotated_genes_right()))
					
					locations = {}
					for location in fusion.locations:
						if(not locations.has_key(location["dataset"])):
							locations[location["dataset"]] = []
						locations[location["dataset"]].append("uid."+location["id"]+"-"+location["left"][0]+":"+str(location["left"][1])+"-"+location["right"][0]+":"+str(location["right"][1]))
					
					for dataset in datasets:
						fh.write("\t")
						if(locations.has_key(dataset)):
							fh.write(",".join(locations[dataset]))
					fh.write("\n")
			fh.close()
	
	def export3(self,filename_prefix,glue=" & "):
		dataset_names = []
		dataset_names_index = {}
		dataset_names_lengths = []
		dataset_names_sized = []
		
		i = 1
		for ds in self.datasets:
			dataset_names.append(ds.name)
			dataset_names_index[ds.name] = i
			dataset_names_lengths.append(len(ds.name))
			dataset_names_sized.append(ds.name+" ("+str(ds.count_fusions())+")")
			i += 1
		dataset_names_length = (len(glue) * (len(self.datasets)-2)) + sum(sorted(dataset_names_lengths)[::-1][:-1])
		
		filename = filename_prefix+"overlap_summary.txt"
		print "Putting output into: "+filename
		fh = open(filename,"w")
		
		for r in range(1,len(self.datasets)):
			#filename = filename_prefix+"matches_"+str(r+1)+".txt"
			#fh = open(filename,"w")
			
			line = "\t"+"\t".join(dataset_names_sized)
			fh.write(line+"\n")
			for x in itertools.combinations(dataset_names,r):
				l = list(x)
				prefix = " & ".join(x)
				spacer = ""#" " * (dataset_names_length - len(prefix))
				
				key = []
				for item in l:
					key.append(str(dataset_names_index[item]))
				key = ".".join(key)
				
				#line_i = "| "+prefix+spacer+" ("+str(self.matrix_tmp[key].count_fusions())+")"
				line_i = prefix+spacer+" ("+str(self.matrix_tmp[key].count_fusions())+")"
				
				i = 0
				for ds in dataset_names:
					line_i += "\t"#" | "
					if(ds in l):
						pass
						#line_i += " " * len(dataset_names_sized[i])
					else:
						key2 = [i+1]
						for item in l:
							key2.append(dataset_names_index[item])
						
						key2 = sorted(key2)
						
						for j in range(len(key2)):
							key2[j] = str(key2[j])
						
						key2 = ".".join(key2)
						
						a = self.matrix_tmp[key2].count_fusions()
						b = self.matrix_tmp[key].count_fusions()
						c = self.datasets[i].count_fusions()
						
						if(b > 0):
							perc1_str = str(a)+"/"+str(b)+" ("+str(round(100.0*a/b,2))+"%)"
						else:
							perc1_str = "0"
						
						if(c > 0):
							perc2_str = str(a)+"/"+str(c)+" ("+str(round(100.0*a/c,2))+"%)"
						else:
							perc2_str = "0"
						
						perc_str = perc1_str + " : " + perc2_str
						line_i += perc_str #" " * (len(dataset_names_sized[i]) - len(perc_str)) + perc_str
					i += 1
				fh.write(line_i+"\n")
			fh.write("\n\n")
		fh.close()
		
		"""
		Create all tables:
		
		_1.txt:
		             | 
		prediction_a | 123
		prediction_b |  22
		prediction_c | 140
		
		_2.txt:
		
		                   | prediction_a (123) | prediction_b (22) | prediction_c (140)
		prediction_a (123) |                    |     123/22 (xxx%) |     123/140 (xxx%)
		prediction_b  (22) |      22/123 (xxx%) |                   |       22/140 (xx%)
		prediction_c (140) |     140/123 (xxx%) |     140/22 (xxx%) |                   
		
		_3.txt
		                            | prediction_a (123) | prediction_b (22) | prediction_c (140)
		prediction_a & prediction_b |                                        |     123/140 (xxx%)
		prediction_a & prediction_c |                    |      22/22 (100%) |                   
		prediction_b & prediction_c |     140/123 (xxx%) |                                       
		
		
		
		
		
		
		4 case:
		
		2:
			| a | b | c | d
		-------------------
		a         +   +   +
		b     +       +   +
		c     +   +       +
		d     +   +   +    
		
		
		3:
				| a | b | c | d
		-------------------
		a & b             +   +
		a & c         +       +
		a & d         +   +
		b & c     +           +
		b & d     +       +
		c & d     +   +
		
		4:
					| a | b | c | d
		-------------------
		a & b & c	| 
		a & b & d	| 
		a & c & d	| 
		b & c & d	| 
		
		
		
					| a | b | c | d
		-------------------
		a & b & c &d| 
		
		"""
	
	"""
	def export_distances(self,filename_prefix="",join="_vs._",suffix=".txt"):
		if(self.distances != False):
			for i in range(len(self.distances)):
				for j in range(len(self.distances[i])):
					filename = filename_prefix+self.datasets[i].name+join+os.path.basename(self.datasets[j].name)+suffix
					print "exporting: "+os.path.basename(filename)
					fh = open(filename,"w")
					fh.write("left_fusion["+self.datasets[i].name+"]\t")
					fh.write("right_fusion["+self.datasets[j].name+"]\t")
					fh.write("eucledian distance\t")
					fh.write("city block distance\n")
					
					distances = self.distances[i][j]
					# sort by eucledian distance
					
					keys = {}
					for distance in distances:
						#print distance
						
						
						if(not keys.has_key(distance["eucledian distance"])):
							keys[distance["eucledian distance"]] = []
						keys[distance["eucledian distance"]].append(distance)
					
					for key in sorted(keys.keys()):
						for distance in keys[key]:
							fh.write(str(distance["fusion_1"].get_left_chromosome())+":")
							fh.write(str(distance["fusion_1"].get_left_break_position())+"-")
							fh.write(str(distance["fusion_1"].get_right_chromosome())+":")
							fh.write(str(distance["fusion_1"].get_right_break_position())+"\t")
							
							fh.write(str(distance["fusion_2"].get_left_chromosome())+":")
							fh.write(str(distance["fusion_2"].get_left_break_position())+"-")
							fh.write(str(distance["fusion_2"].get_right_chromosome())+":")
							fh.write(str(distance["fusion_2"].get_right_break_position())+"\t")
							
							fh.write(str(distance["eucledian distance"])+"\t")
							fh.write(str(distance["city block distance"])+"\n")
					
					fh.close()
	"""


"""
Figure out the amount of comparisons:

with n experiments, we want to overlay k ( 1 < k <= n) experiments in incremental order (usefull for DP):

say n = 4, k = {2,3,4}

because the order does not matter (comparing a with b should give the same result as b with a), we want to find the number of combinations:
we find for k anumber of j = (n! / (n-k)!) comparisons

k = 2:     4! / (4-2)!
"""

