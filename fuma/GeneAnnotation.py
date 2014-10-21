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

import HTSeq

class GeneAnnotation:
	def __init__(self,name):
		self.n = 0
		self.name = name
		self.gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	
	def add_annotation(self,gene,chromosome,start,stop):
		self.n += 1
		self.gas[HTSeq.GenomicInterval(chromosome,start,stop)] += gene
	
	def get_annotations(self,chromosome,position):
		#unique_genes = list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in r])) << weird list construction - only neccesairy using the steps() function
		for annotation in self.gas[HTSeq.GenomicPosition(chromosome,position)]:
			yield annotation
	
	def show_me(self):
		for chromosome_name,chromosome_obj in self.gas.chrom_vectors.items():
			print "Chromosome: "+str(chromosome_name)
			genes = list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in self.gas[HTSeq.GenomicInterval(chromosome_name,0,chromosome_obj['.'].iv.end)].steps()]))
			for gene in genes:
				print " - "+str(gene)
	
	def __len__(self):
		return self.n
	
	def __iter__(self):
		for chromosome_name,chromosome_obj in self.gas.chrom_vectors.items():
			for gene in list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in self.gas[HTSeq.GenomicInterval(chromosome_name,0,chromosome_obj['.'].iv.end)].steps()])):
				yield gene