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
from .FusionDetectionExperiment import FusionDetectionExperiment
from .MergedFusion import MergedFusion

from .Fusion import AD_DIRECTION_REVERSE
from .Fusion import AD_DIRECTION_FORWARD

from .Fusion import STRAND_FORWARD
from .Fusion import STRAND_REVERSE


import os.path,sys,itertools


class ComparisonTriangle:
	logger = logging.getLogger("FuMa::ComparisonTriangle")
	
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
		fh = self.export_list_header()
		
		export_fusions = []# Fusions to be exported after current iteration
		merged_fusions = []# MergedFusions to be used for next iteration
		
		for tmp,fusion in self:
			export_fusions.append(fusion)
		
		n_total = int(round(0.5 * (len(export_fusions) * (len(export_fusions) + 1)) ))
		passed = 0
		previous_percentage = -100.0
		
		self.logger.info("Starting "+str(n_total)+" comparisons for k=1")
		
		for y,fusion_y in self:
			for x,fusion_x in self:
				if y >= x:
					n_total, passed, previous_percentage = self.log_progress(n_total, passed, previous_percentage)
					
					# If they do not belong to the same dataset - i.e. no duplication removal - and if they are the same MergedFusion gene
					if fusion_y and fusion_y.dataset_name not in [tmp['dataset'] for tmp in fusion_x.locations()]:# and fusion_y != fusion_x
						comparison = self.match_fusions(fusion_y, fusion_x)
						
						if comparison != False:
							# Keep is not important - hiding is only useful for for exporting..
							export_fusions[x] = None
							export_fusions[y] = None
							
							merged_fusions.append(comparison)
					passed += 1
		
		n_total, passed, previous_percentage = self.log_progress(n_total, passed, previous_percentage)
		
		self.export_list_chunked(fh,export_fusions)
		
		#@todo put this in some kind of while loop - and add recursion limit to be better safe than sorry..
		while len(merged_fusions) > 0:
			merged_fusions = self.overlay_fusions_recursive(fh,merged_fusions)
		
		if self.args.output != "-":
			fh.close()
	
	def overlay_fusions_recursive(self,fh,merged_fusions):
		n_total = self.num_fusions() * len(merged_fusions)
		passed = 0
		previous_percentage = -100.0
		
		k = len(merged_fusions[0].fusions)
		if k > self.num_fusions():
			raise Exception("Out of bound, reasonable recursion depth has exceeded")
		
		self.logger.info("Starting "+str(n_total)+" comparisons for k="+str(k))
		
		export_fusions = [tmp for tmp in merged_fusions]
		merged_fusions_new = []
		
		for x in range(len(merged_fusions)):
			for y,fusion_y in self:
				merged_fusion_x = merged_fusions[x]
				n_total, passed, previous_percentage = self.log_progress(n_total, passed, previous_percentage)
				
				# - if fusion_y not in merged_fusion ? // if fusion_y.dataset_name not in [tmp['dataset'] for tmp in fusion_x.locations()]
				if fusion_y not in merged_fusion_x.fusions:
					comparison = self.match_fusions(fusion_y, merged_fusion_x)
					
					if comparison != False:
						export_fusions[x] = None
						merged_fusions_new.append(comparison)
						break
				
				passed += 1
		n_total, passed, previous_percentage = self.log_progress(n_total, passed, previous_percentage)
		
		merged_fusions_new = self.prune_duplicates(merged_fusions_new)
		self.export_list_chunked(fh,export_fusions)
		
		return merged_fusions_new
	
	def prune_duplicates(self,merged_fusions):
		"""
		Remove MergedFusion instances with identical Fusion objects
		"""
		
		for y in range(len(merged_fusions)):
			for x in range(y+1,len(merged_fusions)):
				if merged_fusions[x] != None and merged_fusions[y] != None and merged_fusions[x].fusions == merged_fusions[y].fusions:
					# Do some garbage removal?
					merged_fusions[x] = None
		
		return [tmp for tmp in merged_fusions if tmp != None]
	
	def log_progress(self,n_total, passed, previous_percentage):
		# Print percentage - doesn't entirely fit yet
		try:
			percentage = 100.0 * (float(passed) / float(n_total))
		except ZeroDivisionError:
			percentage = 100.0
		if percentage >= previous_percentage + 5.0 or passed == n_total:# Repport each 5%
			self.logger.debug(str(round(percentage,1))+"% completed")
			previous_percentage = percentage
		return n_total, passed, previous_percentage
	
	def export_list_fg(self,fusion,fh):
		if(self.args.acceptor_donor_order_specific_matching and fusion.acceptor_donor_direction == AD_DIRECTION_REVERSE):
			## A-B should be reported as B-A; chr1:123\tchr1:456 as chr1:456-chr1:123
			fh.write(":".join(sorted(list(set([str(gene) for gene in fusion.get_annotated_genes_right2()])))) + "\t")
			fh.write(":".join(sorted(list(set([str(gene) for gene in fusion.get_annotated_genes_left2()])))))
			
			if fusion.spans_a_large_gene():
				fh.write("\tTRUE")
			else:
				fh.write("\tFALSE")
			
			for dataset in self.experiments:
				fh.write("\t")
				strdata = []
				
				for location in fusion.locations():
					if location['dataset'] == dataset.name:
						strdata.append(str(location['id'])+"=chr"+location['right'][0]+':'+str(location['right'][1])+'-chr'+location['left'][0]+':'+str(location['left'][1]))
				
				fh.write(",".join(sorted(strdata)))
		else:
			fh.write(":".join(sorted(list(set([str(gene) for gene in fusion.get_annotated_genes_left2()])))) + "\t")
			fh.write(":".join(sorted(list(set([str(gene) for gene in fusion.get_annotated_genes_right2()])))))
			
			if fusion.spans_a_large_gene():
				fh.write("\tTRUE")
			else:
				fh.write("\tFALSE")
			
			for dataset in self.experiments:
				fh.write("\t")
				strdata = []
				
				for location in fusion.locations():
					if location['dataset'] == dataset.name:
						strdata.append(str(location['id'])+"=chr"+location['left'][0]+':'+str(location['left'][1])+'-chr'+location['right'][0]+':'+str(location['right'][1]))
				
				fh.write(",".join(sorted(strdata)))
			fh.write("\n")
	
	def export_list_header(self):
		if self.args.output == "-":
			fh = sys.stdout
		else:
			fh = open(self.args.output,"w")
		
		fh.write("Left-genes\tRight-genes\t")
		if self.args.long_gene_size > 0:
			fh.write("Spans large gene (>"+str(self.args.long_gene_size)+"bp)")
		else:
			fh.write("Spans large gene (feature disabled)")
		
		for experiment in self.experiments:
			fh.write("\t"+experiment.name)
		fh.write("\n")
		
		return fh
	
	def export_list_chunked(self,fh,chunk_fusions):
		for fusion in chunk_fusions:
			if fusion not in [None, False]:# False means marked as duplicate earlier on
				self.export_list_fg(fusion,fh)
				
				if not isinstance(fusion, Fusion):
					del(fusion)# MergedFusion is not part of any original object, only used for determining the overlap. Get rid of object as soon as possible
	
	def __len__(self):
		return len(self.experiments)
	
	def num_fusions(self):
		n = 0
		for experiment in self.experiments:
			n += len(experiment)
		return n
	
	def __iter__(self):
		i = 0
		for experiment in self.experiments:
			for fusion in experiment:
				yield i,fusion
				i += 1
	
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
		if	self.match_fusion_gene_strands(fusion_1,fusion_2) and \
			self.match_acceptor_donor_direction(fusion_1,fusion_2) and \
			fusion_1.has_annotated_genes() and \
			fusion_2.has_annotated_genes():
			
			# Compare the fusion genes based on their gene names
			if(self.args.matching_method == 'overlap'):
				matches_left  = self.match_overlap(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
				matches_right = self.match_overlap(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
			elif(self.args.matching_method == 'egm'):
				matches_left  = self.match_egm(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
				matches_right = self.match_egm(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
			else:
				matches_left  = self.match_sets(fusion_1.get_annotated_genes_left2(), fusion_2.get_annotated_genes_left2())
				matches_right = self.match_sets(fusion_1.get_annotated_genes_right2(), fusion_2.get_annotated_genes_right2())
			
			if matches_left and matches_right:
				# Fusion only merges with MergedFusion
				#if isinstance(fusion_1, MergedFusion) and isinstance(fusion_2, MergedFusion):
				#	raise Exception("If (A & B) == (C & D), (A & B & C) should have matched before..")
				#	#merged_fusion = fusion_1
				#	#merged_fusion.merge(fusion_2)
				#	#replace_merged_fusions = fusion_2
				#
				# And the  first object is always a Fusion, the second possibly a MergedFusion
				#elif isinstance(fusion_1, MergedFusion) and isinstance(fusion_2, Fusion):
				#	merged_fusion = fusion_1
				#	merged_fusion.add_fusion(fusion_2)
				#
				# And the following can be done cleaner
				#elif isinstance(fusion_1, Fusion) and isinstance(fusion_2, MergedFusion):
				#	merged_fusion = fusion_2
				#	merged_fusion.add_fusion(fusion_1)
				#elif isinstance(fusion_1, Fusion) and isinstance(fusion_2, Fusion):
				#	merged_fusion = MergedFusion()
				#	merged_fusion.add_fusion(fusion_1)
				#	merged_fusion.add_fusion(fusion_2)
				
				if isinstance(fusion_1, Fusion):
					if isinstance(fusion_2, MergedFusion):
						merged_fusion = fusion_2
						merged_fusion.add_fusion(fusion_1)
					elif isinstance(fusion_2, Fusion):
						merged_fusion = MergedFusion()
						merged_fusion.add_fusion(fusion_1)
						merged_fusion.add_fusion(fusion_2)
					else:
						raise Exception("Something went wrong with the object types")
				else:
					raise Exception("Something went wrong with the object types")
				
				# This has to be pre-cached and can not be determined on the fly by a functions,
				# because it requires the type of matching. If you would allow for functions, you could 
				# end up with overlap and egm and subset based matching mixed up.
				merged_fusion.annotated_genes_left = matches_left
				merged_fusion.annotated_genes_right = matches_right
				
				return merged_fusion
		return False
	
	def match_fusion_gene_strands(self,fusion_1,fusion_2):
		if not self.args.strand_specific_matching:
			return True
		else:
			if fusion_1.get_left_strand() == None or fusion_1.get_right_strand() == None or fusion_2.get_left_strand() == None or fusion_2.get_right_strand() == None:
				raise Exception("A fusion gene without an annotated strand was used for strand-specific-matching.\n\n"+fusion_1.__str__()+"\n"+fusion_2.__str__())
			else:
				return fusion_1.get_left_strand() == fusion_2.get_left_strand() and fusion_1.get_right_strand() == fusion_2.get_right_strand()
	
	def match_acceptor_donor_direction(self,fusion_1,fusion_2):
		if(not self.args.acceptor_donor_order_specific_matching):
			return True
		elif(fusion_1.acceptor_donor_direction == None or fusion_2.acceptor_donor_direction == None):
			raise Exception("A fusion gene without an annotated acceptor-donor direction was used for acceptor-donor-order-specific-matching.\n\n"+fusion_1.__str__()+"\n"+fusion_2.__str__())
		else:
			return (fusion_1.acceptor_donor_direction == fusion_2.acceptor_donor_direction)
	
	def match_overlap(self,set1,set2):									#https://docs.python.org/2/library/sets.html
		set1_s = set([str(gene) for gene in set1])
		set2_s = set([str(gene) for gene in set2])
		
		overlap = set1_s.intersection(set2_s)
		if(not overlap):
			return None
		else:
			gene_list = []
			for gene in set1:
				if str(gene) in overlap:
					gene_list.append(gene)
			return gene_list
	
	def match_egm(self,set1,set2):
		set1_s = set([str(gene) for gene in set1])
		set2_s = set([str(gene) for gene in set2])
		if set1_s == set2_s:
			return set1
		else:
			return None
	
	def match_sets(self,superset,subset):								#https://docs.python.org/2/library/sets.html
		subset_s = set([str(gene) for gene in subset])
		superset_s = set([str(gene) for gene in superset])
		
		if(len(subset_s) > len(superset_s)):
			return self.match_sets(subset,superset)						# Gene names have to be provided as sets
		elif(subset_s.issubset(superset_s)):
			return subset
		else:
			return None
