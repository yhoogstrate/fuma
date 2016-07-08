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

"""
Example of the BED format:

#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
chr1	67075869	67163158	NM_207014	0	-	67075923	67163102	0	10	198,203,195,156,140,157,113,185,175,226,	0,2870,9885,24548,33771,37182,53555,55630,67602,87063,
chr1	67051159	67163158	NM_024763	0	-	67052400	67163102	0	17	1292,157,227,99,122,158,152,87,203,195,156,140,157,113,185,175,226,	0,9472,13931,14923,20696,21102,22737,24821,27580,34595,49258,58481,61892,78265,80340,92312,111773,
chr1	8335050	8800286	NM_001042681	0	-	8337733	8638943	0	23	2717,181,147,721,223,1379,114,162,200,93,163,81,99,100,125,49,105,97,106,126,71,469,481,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,113521,142659,145001,156222,188809,204070,205013,262156,271905,303568,464755,
chr1	8335050	8800286	NM_012102	0	-	8337733	8638943	0	24	2717,181,147,721,223,1379,114,162,200,93,163,81,99,100,125,49,105,97,106,126,71,469,185,481,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,113521,142659,145001,156222,188809,204070,205013,262156,271905,303568,439999,464755,
chr1	8335050	8406334	NM_001042682	0	-	8337733	8346780	0	13	2717,181,147,721,223,1379,114,162,200,93,163,81,127,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,71157,
"""

import logging
import sys

from .Gene import Gene
from .GeneAnnotation import GeneAnnotation

import HTSeq

class ParseBED(GeneAnnotation):
	logger = logging.getLogger("FuMa::ParseBED")
	
	def __init__(self,filename,name,long_gene_size):
		GeneAnnotation.__init__(self,name)
		self.long_gene_size = long_gene_size
		self.parse(filename)
	
	def parse(self,filename):
		self.logger.info('Parsing BED file: '+str(filename))
		
		with open(filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
		
		self.logger.debug('Size of Gene Annotation: '+str(len(self)))
	
	def cleanup_chr_name(self,chr_name):
		"""Given the large number of fusion genes, we remove all 'chr'
		prefixes because they add 6 bytes per fusion gene. They can be
		added again using the getter functions.
		"""
		
		chr_name = chr_name.strip()
		return chr_name[3:] if chr_name[0:3] == "chr" else chr_name
	
	def parse_line(self,line):
		line = line.split("\t")
		if(len(line) >= 4):
			
			start = int(line[1])
			stop = int(line[2])
			length = abs(stop - start)
			is_long_gene = (length > self.long_gene_size) and (self.long_gene_size != 0)
			self.add_annotation(Gene(line[3], is_long_gene), self.cleanup_chr_name(line[0]), start, stop )
		
		#@ deprecated - exon-type BED files
		#self.index[line[3]].append([line[0],int(line[1]),int(line[2])])
	
	"""
	@deprecated functions
	@ - conversion of exon-type BED to gene-type BED will suffer from duplicated genes on the same chromosome
	
	def merge_exons_to_gene(self,gene_id,exons):
		pos_min = exons[0][1]
		pos_max = exons[0][2]
		
		chromosomes = []
		
		for exon in exons[1:]:
			if(exon[0] != exons[0][0]):
				chromosomes.append(exon[0])
			else:
				if(exon[1] < pos_min):
					pos_min = exon[1]
				if(exon[2] > pos_max):
					pos_max = exon[2]
		
		if(len(chromosomes) > 0):
			self.logger.error("Gene '"+gene_id+"' ("+exons[0][0] + ") has multiple entries AND is also located at other chromosomes: '" + "' & '".join(list(set(chromosomes))) + "'. Please double check your BED file.")
		
		return [exons[0][0],pos_min,pos_max]
	
	def push_index(self):
		for gene_id in self.index.keys():
			items = self.index[gene_id]
			if(len(items) > 1):
				gene = self.merge_exons_to_gene(gene_id,items)
				#self.logger.warning("Gene '"+gene_id+"' has "+str(len(items))+" annotations instead of 1, and is therefore merged into: "+gene[0]+":"+str(gene[1])+"-"+str(gene[2]))
			else:
				gene = items[0]
			
			self.add_annotation(Gene(gene_id),gene[0],int(gene[1]),int(gene[2]))
	"""
