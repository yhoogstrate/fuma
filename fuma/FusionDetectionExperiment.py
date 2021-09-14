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

import logging,sys,fuma,datetime

from fuma import Fusion
from fuma.Fusion import AD_DIRECTION_REVERSE
from fuma.Fusion import AD_DIRECTION_FORWARD

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE

class FusionDetectionExperiment:
	logger = logging.getLogger("FuMa::FusionDetectionExperiment")
	
	def __init__(self,name):
		self.name = name
		
		self.genes_spanning_left_junction = None
		self.genes_spanning_right_junction = None
		
		self.flush()
	
	def __del__(self):
		for fusion in self:
			fusion.prepare_deletion()
	
	def add_fusion(self,fusion):
		# Add left location
		left_chr = fusion.get_left_chromosome(False)
		left_pos = fusion.get_left_break_position()
		
		# Add right location
		right_chr = fusion.get_right_chromosome(False)
		right_pos = fusion.get_right_break_position()
		
		###################### new type of indexing ####################
		## ensure that chr_left < chr_right
		if(left_chr not in self.index):
			self.index[left_chr] = {}
		
		if(right_chr not in self.index[left_chr]):
			self.index[left_chr][right_chr] = []
		
		self.index[left_chr][right_chr].append(fusion)
		################################################################
		
		self.n += 1
	
	def show_me(self):
		print(self.__str__())
	
	def __str__(self):
		out = "---------------------\n"
		out += "Showing experiment: "+self.name+"\n\n"
		for fusion in self.__iter__():
			if(fusion):# Duplicates are flagged as False/None
				out += fusion.__str__()+"\n"
		out += "---------------------\n"
		
		return out
	
	def export_to_CG_Junctions_file(self,filename):
		if(filename == "-"):
			fh = sys.stdout
		else:
			fh = open(filename,"w")
		
		fh.write("#ASSEMBLY_ID	???\n")
		fh.write("#SOFTWARE_VERSION	FuMa v"+fuma.__version__+"\n")
		fh.write("#GENERATED_BY	FuMa\n")
		fh.write("#GENERATED_AT	"+str(datetime.datetime.utcnow())+"\n")
		fh.write("#FORMAT_VERSION	2\n")
		fh.write("#GENOME_REFERENCE	???	build	??\n")
		fh.write("#SAMPLE	???\n")
		fh.write("#TYPE	JUNCTIONS\n")
		fh.write("#DBSNP_BUILD	dbSNP	build	???\n")
		fh.write("#GENE_ANNOTATIONS	???	build	???\n")
		fh.write("\n")
		fh.write(">Id	LeftChr	LeftPosition	LeftStrand	LeftLength	RightChr	RightPosition	RightStrand	RightLength	StrandConsistent	Interchromosomal	Distance	DiscordantMatePairAlignments	JunctionSequenceResolved	TransitionSequence	TransitionLength	LeftRepeatClassification	RightRepeatClassification	LeftGenes	RightGenes	XRef	DeletedTransposableElement	KnownUnderrepresentedRepeat	FrequencyInBaselineGenomeSet	AssembledSequence	EventId	Type	RelatedJunctions\n")
		
		fid = 1
		
		for fusion in self.__iter__():
			if(fusion != False):# Duplicates are flagged as False
				fh.write(str(fid)+"	")
				
				fh.write(fusion.get_left_chromosome()+"	")
				fh.write(str(fusion.get_left_break_position())+"	")
				
				if(fusion.left_strand == STRAND_FORWARD):
					fh.write('+')
				elif(fusion.left_strand == STRAND_REVERSE):
					fh.write('-')

				fh.write("	101	")
				fh.write(fusion.get_right_chromosome()+"	")
				fh.write(str(fusion.get_right_break_position())+"	")
				
				if(fusion.left_strand == STRAND_FORWARD):
					fh.write('+')
				elif(fusion.left_strand == STRAND_REVERSE):
					fh.write('-')
				
				fh.write("	101	")
				strand_consistent = (fusion.left_strand == fusion.right_strand)
				interchromosomal = fusion.is_interchromosomal()#(fusion.get_left_chromosome() != fusion.get_right_chromosome())
				distance = str(fusion.get_distance())
				
				if(strand_consistent):
					fh.write("Y	")
				else:
					fh.write("N	")
				
				if(interchromosomal):
					fh.write("Y	")
				else:
					fh.write("N	")
				
				if(distance != "-1"):
					fh.write(distance)
				fh.write("\t")
				fh.write("20	Y	")
				
				# Sequences are not being stored as it consumes too much memory
				fh.write("\t")
				
				fh.write("			")
				
				fh.write(":".join(fusion.get_annotated_genes_left(True).keys())+"	")
				fh.write(":".join(fusion.get_annotated_genes_right(True).keys())+"	")
				
				fh.write("			1.0		"+str(fid)+"	complex	"+str(fusion.locations)+"\n")
				
				fid += 1
		
		if(filename != "-"):
			fh.close()
	
	def export_to_list(self,fh,order,blacklist,args):
		"""
		Exports to a tabular file of the following syntax:
		Left-parner(s) \t right-parner(s) \t dection-method 1     \t detection-method 2
		TMPRSS2        \t ERG             \t chr21:...-chrr21:... \t chr21:...-chr21
		"""
		for fusion in self.__iter__():
			cur_datasets = fusion.dataset_name.split('_vs._')
			
			check = True
			for initial_fusion in fusion.matches:
				if(initial_fusion in blacklist):
					check = False
					break
			
			# Duplicates are flagged as False
			if(fusion != False and fusion.get_dataset_statistics()[1] == 0 and check):
				if(fusion.acceptor_donor_direction == AD_DIRECTION_REVERSE and args.acceptor_donor_order_specific_matching):
					# A-B should be reported as B-A; chr1:123\tchr1:456 as chr1:456-chr1:123
					
					fh.write(":".join(sorted(fusion.get_annotated_genes_right(True).keys()))+"	")
					fh.write(":".join(sorted(fusion.get_annotated_genes_left(True).keys()))+"")
					
					if fusion.spans_a_large_gene():
						fh.write("\tTRUE")
					else:
						fh.write("\tFALSE")
					
					for dataset in order:
						strdata = []
						try:
							i = cur_datasets.index(dataset)
						except:
							i = -1
						
						fh.write("\t")
						if(i > -1):
							for loc in fusion.locations():
								if(loc['dataset'] == dataset):
									strdata.append(str(loc['id'])+"=chr"+loc['right'][0]+':'+str(loc['right'][1])+'-chr'+loc['left'][0]+':'+str(loc['left'][1]))
							fh.write(",".join(sorted(strdata)))
				else:
					fh.write(":".join(sorted(fusion.get_annotated_genes_left(True).keys()))+"	")
					fh.write(":".join(sorted(fusion.get_annotated_genes_right(True).keys()))+"")
					
					if fusion.spans_a_large_gene():
						fh.write("\tTRUE")
					else:
						fh.write("\tFALSE")
					
					for dataset in order:
						strdata = []
						try:
							i = cur_datasets.index(dataset)
						except:
							i = -1
						
						fh.write("\t")
						if(i > -1):
							for loc in fusion.locations():
								if(loc['dataset'] == dataset):
									strdata.append(str(loc['id'])+"=chr"+loc['left'][0]+':'+str(loc['left'][1])+'-chr'+loc['right'][0]+':'+str(loc['right'][1]))
							fh.write(",".join(sorted(strdata)))
				fh.write("\n")
	
	def annotate_genes(self,gene_annotation):
		self.annotate_genes_left(gene_annotation)
		self.annotate_genes_right(gene_annotation)
	
	def annotate_genes_left(self,gene_annotation):
		if(not self.genes_spanning_left_junction):
			self.logger.debug("Annotating genes on the left junction: "+self.name+" - "+gene_annotation.name)
			
			for fusion in self.__iter__():
				if(fusion.annotated_genes_left == None):				# if object is not set, make it an empty list
					fusion.annotated_genes_left = []
				
				for gene in gene_annotation.get_annotations(fusion.get_left_chromosome(),fusion.get_left_break_position()):
					fusion.annotated_genes_left.append(gene)
			
			self.genes_spanning_left_junction = [gene_annotation]
	
	def annotate_genes_right(self,gene_annotation):
		if(not self.genes_spanning_right_junction):
			self.logger.debug("Annotating genes on the right junction: "+self.name+" - "+gene_annotation.name)
			
			for fusion in self:
				if(fusion.annotated_genes_right == None):				# if object is not set, make it an empty list
					fusion.annotated_genes_right = []
				
				for gene in gene_annotation.get_annotations(fusion.get_right_chromosome(),fusion.get_right_break_position()):
					fusion.annotated_genes_right.append(gene)
			
			self.genes_spanning_right_junction = [gene_annotation]
	
	def __iter__(self):
		""" Return all fusions (non-indexed but sorted on chr,chr) as iterator
		"""
		for key1 in sorted(self.index.keys()):
			chromosome_left = self.index[key1]
			
			for key2 in sorted(chromosome_left.keys()):
				chromosome_right = chromosome_left[key2]
				
				for fusion in chromosome_right:
					yield fusion
	
	def __getitem__(self,i):
		"""
		This allows to acces fusions by their id.
		For:
		fusions = FusionDetectionExperiment
		
		we can obtain the fourth fusion by:
		
		fusion_4th = fusions[3]# Use 0 based counting
		"""
		j = 0
		for fusion in self:
			if(i == j):
				return fusion
			else:
				j += 1
		
		return None
	
	def remove_duplicates(self,args):
		"""
		- First create a table of those that overlap
		- Then create merged entries based on the overlap matrix
		"""
		if(not self.genes_spanning_left_junction or not self.genes_spanning_right_junction):
			raise Exception("Gene annotations on dataset '"+self.name+"' were not found")
		else:
			old_count = len(self)
			if(self.name.find("vs.") == -1):
				self.logger.info("Duplication removal: "+self.name+" ("+str(old_count)+" fusions)")
		
		unique_fusions = []
		
		if(args.matching_method in ["overlap","subset","egm"]):
			from .CompareFusionsBySpanningGenes import CompareFusionsBySpanningGenes
			overlap = CompareFusionsBySpanningGenes(False,False,args)
		else:
			raise Exception("Unknown overlap method for removing duplicates: '"+args.matching_method+"' for dataset "+self.name)
		
		stats_duplicates = 0
		stats_non_gene_spanning = 0
		
		fusions_to_add = []
		
		for chromosome_left in self.index.items():
			for chromosome_right in chromosome_left[1].items():
				
				all_fusions = chromosome_right[1]
				n = len(all_fusions)
				
				queue = range(n)
				while(len(queue) > 0):
					duplicates = []
					for i in queue:
						fusion_1 = all_fusions[i]
						if(fusion_1):
							is_duplicate = False
							if(len(fusion_1.get_annotated_genes_left(False)) == 0 or len(fusion_1.get_annotated_genes_right(False)) == 0):
								stats_non_gene_spanning += 1
								all_fusions[i] = False
							else:
								for j in range(i+1,n):
									fusion_2 = all_fusions[j]
									if(fusion_2):
										match = overlap.match_fusions(fusion_1,fusion_2,False)
										
										if(match):
											merged_matches = fusion_1.matches | fusion_2.matches
											
											fusion_1.matches = merged_matches
											fusion_1.acceptor_donor_direction = match.acceptor_donor_direction
											fusion_1.left_strand = match.left_strand
											fusion_1.right_strand = match.right_strand
											fusion_1.annotated_genes_left  = match.annotated_genes_left
											fusion_1.annotated_genes_right  = match.annotated_genes_right
											
											all_fusions[i] = fusion_1
											all_fusions[j] = False
											is_duplicate = True
											
											match.prepare_deletion()
											del(match)
								
								if(is_duplicate):
									duplicates.append(i)
								else:
									unique_fusions.append(fusion_1)
					queue = duplicates
				
				for fusion in all_fusions:
					if(fusion):
						fusions_to_add.append(fusion)
		
		self.flush()
		for fusion in fusions_to_add:
			self.add_fusion(fusion)
		
		if(self.name.find("vs.") == -1):
			self.logger.debug("* Full: "+str(old_count))
			self.logger.debug("* Gene-spanning: "+str(old_count-stats_non_gene_spanning))
			self.logger.debug("* Unique: "+str(len(self)))
		
		return len(self)
	
	def __len__(self):
		return self.n
	
	def flush(self):
		self.n = 0
		
		self.n_matches_exp_1 = None
		self.n_matches_exp_2 = None
		
		self.index = {}
