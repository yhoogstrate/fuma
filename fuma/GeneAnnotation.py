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

#import gffutils
import HTSeq
import logging

class GeneAnnotation:
	"""Gene annotation is a virtual reference genome. It's only being
	used the map Genes to in order to index them quickly. 
	"""
	logger = logging.getLogger("FuMa::GeneAnnotation")
	
	def __init__(self,name):
		self.n = 0
		self.name = name
		
		# list(db.region(region=('2L', 9277, 10000), completely_within=True))
		#self.gas2 = gffutils.create_db(gtf, dbfn=db_file)
		self.gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
	
	def add_annotation(self,gene,chromosome,start,stop):
		#self.logger.debug("Adding annotation "+str(self.n)+": "+chromosome+":"+str(start)+"-"+str(stop)+" = "+str(gene))
		self.gas[HTSeq.GenomicInterval(chromosome,start,stop)] += gene
		self.n += 1
	
	def get_annotations(self,chromosome,position):
		#unique_genes = list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in r])) << weird list construction - only neccesairy using the steps() function
		for annotation in self.gas[HTSeq.GenomicPosition(chromosome,position)]:
			yield annotation
	
	def __str__(self):
		out = "[ Gene annotation: "+str(self.name)+" (genes: "+str(len(self))+")]"
		for chromosome_name,chromosome_obj in self.gas.chrom_vectors.items():
			out += "Chromosome: "+str(chromosome_name)+"\n"
			genes = list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in self.gas[HTSeq.GenomicInterval(chromosome_name,0,chromosome_obj['.'].iv.end)].steps()]))
			for gene in genes:
				out += " - "+str(gene)+"\n"
		
		return out
	
	def __len__(self):
		return self.n
	
	def __iter__(self):
		for chromosome_name,chromosome_obj in self.gas.chrom_vectors.items():
			for gene in list(reduce(lambda s1, s2: s1 | s2, [x[1] for x in self.gas[HTSeq.GenomicInterval(chromosome_name,0,chromosome_obj['.'].iv.end)].steps()])):
				yield gene
	
	#def show_me(self):
	#	print self.__str__()