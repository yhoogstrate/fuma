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

from .Readers import *

from .ParseBED import ParseBED
from .CompareFusionsBySpanningGenes import CompareFusionsBySpanningGenes


import os.path,sys,itertools


class OverlapComplex:
	logger = logging.getLogger("FuMa::OverlapComplex")
	
	def __init__(self):
		self.datasets = []
		self.dataset_names = []
		self.matches = False
		self.matches_total = {}
	
	def add_experiment(self,experiment):
		self.datasets.append(experiment)
		self.dataset_names.append(experiment.name)
		self.matches_total[str(len(self.datasets))] = len(experiment)
	
	#def create_matrix(self,n):
	#	matrix = []
	#	for dataset in self.datasets:
	#		matrix.append([])
	#	
	#	n = len(self.datasets)-1
	#	while(n >= 1):
	#		for i in range(n):
	#			matrix[n].append(i)
	#		n -= 1
	#	
	#	return matrix
	
	def create_keys(self,combination):
		keys = [combination[:-1],combination[-1:],combination]
		
		for i in range(len(keys)):
			for j in range(len(keys[i])):
				keys[i][j] = str(keys[i][j])
			
			keys[i] = ".".join(keys[i])
		
		return keys
	
	#def overlay_fusions(self,sparse=True,export_dir=False,output_format="list",matching_method="overlap",strand_specific_matching=False):
	## @todo mv sparse=True into sparse
	## @todo mv export_dir=False into export_dir
	## @todo mv args=None into args
	def overlay_fusions(self,sparse=True,export_dir=False,args=None):
		"""
		The SPARSE variable should only be True if the outpot format
		is 'summary', because all the overlap objects are removed.
		This makes the algorithm much more effictent (reduces space
		complexity from 0.5(n^2) => 2n).
		"""
		n = len(self.datasets)
		
		self.logger.info("Determining the overlap of fusion genes in "+str(n)+" datasets")
		
		self.matrix_tmp = {}
		
		for i in range(len(self.datasets)):
			self.matrix_tmp[str(i+1)] = self.datasets[i]
		
		#comparisons = self.find_combination_table(n)
		if(args.format=="list" and export_dir != False):
			if args.long_gene_size > 0:
				large_genes = "Spans large gene (>"+str(args.long_gene_size)+"bp)"
			else:
				large_genes = "Spans large gene (feature disabled)"
			
			export_dir.write("Left-genes\tRight-genes\t"+large_genes+"\t"+"\t".join(self.dataset_names)+"\n")
		
		ri = 0
		for r in self.find_combination_table(len(self.datasets)):
			r_0 = self.find_combination_table_r_i(len(self.datasets),ri,0)
			
			# First cleanup the memory - reduces space complexity from 0.5(n^2) => 2n. In addition, memory should decrease in time
			dont_remove = []
			matches_this_iteration = set([])
			
			#for c in r:
				#keys = self.create_keys(c)
				
				#dont_remove.append(keys[0])
				#dont_remove.append(keys[1])
			
			#if(args.format != "list"):
				#for candidate in self.matrix_tmp.keys():
					#if candidate not in dont_remove:
						#del(self.matrix_tmp[candidate])
			
			# Then run analysis
			for c in r:
				keys = self.create_keys(c)
				
				comparison = CompareFusionsBySpanningGenes(self.matrix_tmp[keys[0]],self.matrix_tmp[keys[1]],args)
				matches = comparison.find_overlap()
				matches_this_iteration = matches_this_iteration | matches[3]
				
				if(not sparse and export_dir):
					if(args.format=="extensive"):
						matches[0].export_to_CG_Junctions_file(export_dir+"/"+matches[0].name+".CG-junctions.txt")
				
				self.matrix_tmp[keys[2]] = matches[0]
				self.matches_total[keys[2]] = len(matches[0])
			
			if(args.format=="list"):# Write those that are not marked to go to the next iteration to a file
				if(len(r_0) > 2):
					for export_key in self.find_combination_table_r(len(self.datasets),ri-1):#previous_comparisons:#comparisons[ri-1]:
						export_key = [str(x) for x in export_key]
						export_key = '.'.join(export_key)
						
						self.matrix_tmp[export_key].export_to_list(export_dir,self.dataset_names,matches_this_iteration,args)
						del(self.matrix_tmp[export_key]) ## if this was once in a list to be removed, remove...
				else:
					for export_key in [str(i+1) for i in range(len(self.datasets))]:
						self.matrix_tmp[export_key].export_to_list(export_dir,self.dataset_names,matches_this_iteration,args)
						#del(self.matrix_tmp[export_key]) ## if this was once in a list to be removed, remove...
			
			ri += 1
		
		if(args.format == "list" and export_dir != False):
			export_key = '.'.join([str(x) for x in r_0])
			self.matrix_tmp[export_key].export_to_list(export_dir,self.dataset_names,set([]),args) ## if this was once in a list to be removed, remove...?
		
		return matches
	
	def find_combination_table(self,n):
		in_list = range(1,n+1)
		
		for r in range(2,len(in_list)+1):
			yield (list(j) for j in itertools.combinations(in_list,r))
	
	def find_combination_table_r(self,n,r):
		in_list = range(1,n+1)
		r += 2
		return (list(j) for j in itertools.combinations(in_list,r))
	
	def find_combination_table_r_i(self,n,r,i):
		in_list = range(1,n+1)
		r += 2
		return list(list(itertools.combinations(in_list,r))[i])
	
	def set_annotation(self,arg_gene_annotation):
		self.gene_annotation = arg_gene_annotation
	
	def export(self,filename_prefix="",join="_vs._",suffix=".txt"):
		"""
		Exports the very final result of comparison as a table (only those fusions present in all datasets)
		"""
		for i in range(len(self.matches)):
			for j in range(len(self.matches[i])):
				filename = filename_prefix+self.datasets[i].name+join+os.path.basename(self.datasets[j].name)+suffix
				fh = open(filename,"w")
				
				fh.write("["+self.datasets[i].name+"]-position\t")
				fh.write("["+self.datasets[i].name+"]-left-junction-associated-genes\t")
				fh.write("["+self.datasets[i].name+"]-right-junction-associated-genes\t")
				fh.write("["+self.datasets[j].name+"]-position\t")
				fh.write("["+self.datasets[j].name+"]-left-junction-associated-genes\t")
				fh.write("["+self.datasets[j].name+"]-right-junction-associated-genes\n")
				
				for match in self.matches[i][j]:
					fh.write(str(match[0].get_left_chromosome(True))+":"+str(match[0].get_left_break_position())+"-"+str(match[0].get_right_chromosome(True))+":"+str(match[0].get_right_break_position())+"\t")
					fh.write(";".join(match[0].get_annotated_genes_left())+"\t")
					fh.write(";".join(match[0].get_annotated_genes_right())+"\t")
					fh.write(str(match[1].get_left_chromosome(True))+":"+str(match[1].get_left_break_position())+"-"+str(match[1].get_right_chromosome(True))+":"+str(match[1].get_right_break_position())+"\t")
					fh.write(";".join(match[1].get_annotated_genes_left())+"\t")
					fh.write(";".join(match[1].get_annotated_genes_right())+"\n")
				fh.close()
	
	def export_summary(self,filename,glue=" & "):
		"""
		Creates all tables (summary dataformat):
		
		table with 1*1 matches
		
		                   | prediction_a (123) | prediction_b (22) | prediction_c (140)
		prediction_a (123) |                    |     123/22 (xxx%) |     123/140 (xxx%)
		prediction_b  (22) |      22/123 (xxx%) |                   |       22/140 (xx%)
		prediction_c (140) |     140/123 (xxx%) |     140/22 (xxx%) |                   
		
		table with 1*2 matches
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
		
		dataset_names = []
		dataset_names_index = {}
		dataset_names_lengths = []
		dataset_names_sized = []
		
		for i in range(0,len(self.datasets)):
			ds = self.datasets[i]
			dataset_names.append(ds.name)
			dataset_names_index[ds.name] = i+1
			dataset_names_sized.append(ds.name+" ("+str(len(ds))+")")
		
		if(filename == "-"):
			fh = sys.stdout
		else:
			fh = open(filename,"w")
		
		for r in range(1,len(self.datasets)):							# r is the number of datasets merged in the left column
			line = "\t"+"\t".join(dataset_names_sized)
			fh.write(line+"\n")
			
			for x in itertools.combinations(dataset_names,r):
				l = list(x)
				
				key_vertical = ".".join([str(dataset_names_index[item]) for item in x])
				line_i = " & ".join(x)+" ("+str(self.matches_total[key_vertical])+")"
				
				for key_horizontal in range(len(dataset_names)):
					ds = dataset_names[key_horizontal]
					
					line_i += "\t"
					
					if(ds not in l):
						key_merged = [key_horizontal+1]+[dataset_names_index[item] for item in l]
						key_merged = ".".join([str(x) for x in sorted(key_merged)])
						
						a = self.matches_total[key_merged]
						b = self.matches_total[key_vertical]
						c = len(self.datasets[key_horizontal])
						
						if(b > 0):
							perc1_str = str(a)+"/"+str(b)+" ("+str(round(100.0*a//b,2))+"%)"
						else:
							perc1_str = "0"
						
						if(c > 0):
							perc2_str = str(a)+"/"+str(c)+" ("+str(round(100.0*a//c,2))+"%)"
						else:
							perc2_str = "0"
						
						line_i += perc1_str + " : " + perc2_str
				fh.write(line_i+"\n")
			fh.write("\n\n")
		
		if(filename != "-"):
			fh.close()
