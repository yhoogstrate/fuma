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
Small example of the BED format:

#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
chr1	67075869	67163158	NM_207014	0	-	67075923	67163102	0	10	198,203,195,156,140,157,113,185,175,226,	0,2870,9885,24548,33771,37182,53555,55630,67602,87063,
chr1	67051159	67163158	NM_024763	0	-	67052400	67163102	0	17	1292,157,227,99,122,158,152,87,203,195,156,140,157,113,185,175,226,	0,9472,13931,14923,20696,21102,22737,24821,27580,34595,49258,58481,61892,78265,80340,92312,111773,
chr1	8335050	8800286	NM_001042681	0	-	8337733	8638943	0	23	2717,181,147,721,223,1379,114,162,200,93,163,81,99,100,125,49,105,97,106,126,71,469,481,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,113521,142659,145001,156222,188809,204070,205013,262156,271905,303568,464755,
chr1	8335050	8800286	NM_012102	0	-	8337733	8638943	0	24	2717,181,147,721,223,1379,114,162,200,93,163,81,99,100,125,49,105,97,106,126,71,469,185,481,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,113521,142659,145001,156222,188809,204070,205013,262156,271905,303568,439999,464755,
chr1	8335050	8406334	NM_001042682	0	-	8337733	8346780	0	13	2717,181,147,721,223,1379,114,162,200,93,163,81,127,	0,3015,3696,5792,7360,7708,9359,10279,11652,12342,13408,70323,71157,
"""

import sys

from Gene import Gene
from GeneAnnotation import GeneAnnotation

import HTSeq

class ParseBED(GeneAnnotation):
	def __init__(self,filename,name):
		GeneAnnotation.__init__(self,name)
		#self.name = name
		#self.annotations = {}
		#self.annotations_left_indexed = {}
		
		self.parse(filename)
	
	def parse(self,filename):
		with open(filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
		
		#self.index()
	
	def parse_line(self,line):
		line = line.split("\t")
		self.add_annotation(Gene(line[3]),line[0],int(line[1]),int(line[2]))
	
	"""
	def add_gene_annotation(self,name,chromosome,start,stop):
		#Please do some appropriate indexing!
		
		#No index
		
		name = name.split(".")[0]
		
		chromosome = chromosome.replace("chr","")
		if(chromosome.isdigit()):
			chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		
		self.annotations[name] = [chromosome,start,stop]
	"""
	
	"""
	def index(self):
		# Smart index:
		## index[chr1] = [1,2,3,3,4] <- sorted on start_position
		
		tmp_index_left = {}
		
		for name in self.annotations.keys():
			annotation = self.annotations[name]
			chromosome = annotation[0]
			start = annotation[1]
			stop = annotation[2]
			
			if(tmp_index_left.has_key(chromosome) == False):
				tmp_index_left[chromosome] = {}
			
			if(tmp_index_left[chromosome].has_key(start) == False):
				tmp_index_left[chromosome][start] = []
			
			tmp_index_left[chromosome][start].append({"name":name,"chr":chromosome,"start":start,"stop":stop})
		
		
		for chromosome in tmp_index_left.keys():
			self.annotations_left_indexed[chromosome] = []
			for start_position in sorted(tmp_index_left[chromosome].keys()):
				for gene_annotation in tmp_index_left[chromosome][start_position]:
					self.annotations_left_indexed[chromosome].append(gene_annotation)
	"""
	
	"""
	def find_overlap2(self,position,gene):
		#print "find overlap between: ",position,gene,":\t",(position >= gene["start"] and position <= gene["stop"])
		return (position >= gene["start"] and position <= gene["stop"])
	
	def convert_left_locations_to_genes(self,dataset):
		print "Error wrong call; outdated software version!"
		
		import sys
		sys.exit()
		
		#Gene sorting:
		#[  ]
		  #[    ]
		   #[ ]    <- pitfall with loops
		    #[ ]   <- pitfall with loops
		    #[  ]
		    #|
		
		
		#1st:
		 #scan all genes with GENE-START <= fusion.pos
		 #remember pos_1
		#2nd:
		 #scan UNTIL gene with GENE-START > fusion.pos  (while GENE-START <= fusion.pos
		
		#next iteration start from pos_1
		
		
		
		
		
		
		#Possibility:
		
		#ref:
		      #[     ]     <- 1st gene
		 #|                <- 1st event
		
		
		#ref:
		      #[     ]     <- 1st gene
		        #|         <- 1st event
		
		#ref:
		      #[     ]     <- 1st gene
		               #|  <- 1st event
		
		if(not dataset.converted_to_genes_left()):
			print " - Converting LEFT breakpoints to GENES: "+dataset.name
			for chromosome in dataset.get_fusions_indexed_left():
				i = 0
				if(self.annotations_left_indexed.has_key(chromosome["name"])):
					genes = self.annotations_left_indexed[chromosome["name"]]
					for fusion in chromosome["fusions"]:
						current = fusion.get_annotated_genes_left()
						
						if(not current):
							fusion.annotate_genes_left([])
						
						for k in range(0,len(genes)-1):
							gene = genes[k]
							
							if self.find_overlap2(fusion.get_left_break_position(),gene):
								current = fusion.get_annotated_genes_left()
								if(gene["name"] not in current):
									current.append(gene["name"])
								else:
									print gene["name"],chromosome["name"]
								fusion.annotate_genes_left(current)
							
							k += 1
						
				
				
				dataset.is_converted_to_genes = True
		
		print "---"
	
	def convert_right_locations_to_genes(self,dataset):
		print "Error wrong call; outdated software version!"
		
		import sys
		sys.exit()
		
		if(not dataset.converted_to_genes_right()):
			print " - Converting RIGHT breakpoints to GENES: "+dataset.name
			for chromosome in dataset.get_fusions_indexed_right():
				i = 0
				
				if(self.annotations_left_indexed.has_key(chromosome["name"])):
					genes = self.annotations_left_indexed[chromosome["name"]]
					
					for position in chromosome["fusions"]:
						for fusion in position:
							current = fusion.get_annotated_genes_right()
							
							if(not current):
								fusion.annotate_genes_right([])
							
							while(i < len(genes) and fusion.get_right_break_position() < genes[i]["start"]):
								i += 1
							
							k = i
							while(k < len(genes) and fusion.get_right_break_position() >= genes[k]["start"]):
								gene = genes[k]
								if self.find_overlap2(fusion.get_right_break_position(),gene):
									current = fusion.get_annotated_genes_right()
									if(not current):
										current = []
									current.append(gene["name"])
									fusion.annotate_genes_right(current)
								
								k += 1
				
				dataset.is_converted_to_genes = True
	"""
