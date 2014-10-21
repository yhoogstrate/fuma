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
format:

#name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	proteinID	alignID
uc001aaa.2	chr1	+	1115	4121	1115	1115	3	1115,2475,3083,	2090,2584,4121,		uc001aaa.2
uc009vip.1	chr1	+	1115	4272	1115	1115	2	1115,2475,	2090,4272,		uc009vip.1
uc001aab.2	chr1	-	4268	14764	4268	4268	10	4268,4832,5658,6469,6716,7095,7468,7777,8130,14600,	4692,4901,5810,6628,6918,7231,7605,7924,8242,14764,		uc001aab.2
uc009viq.1	chr1	-	4268	19221	4268	4268	7	4268,5658,6469,6720,7468,14600,19183,	4692,5810,6628,6918,7924,14754,19221,		uc009viq.1
"""

class ParseUCSCTable:
	def __init__(self,filename,name):
		self.filename = filename
		self.name = name
		self.annotations = {}
		self.annotations_left_indexed = {}
		
		self.parse()
	
	def parse(self):
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
		
		self.index()
	
	def parse_line(self,line):
		line = line.split("\t")
		self.add_gene_annotation(line[0],line[1],line[3],line[4])
	
	def add_gene_annotation(self,name,chromosome,start,stop):
		"""
		Please do some appropriate indexing!
		"""
		
		#No index
		
		name = name.split(".")[0]
		
		chromosome = chromosome.replace("chr","")
		if(chromosome.isdigit()):
			chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		
		self.annotations[name] = [chromosome,start,stop]
	
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
	
	
	def find_overlap2(self,position,gene):
		#print "find overlap between: ",position,gene,":\t",(position >= gene["start"] and position <= gene["stop"])
		return (position >= gene["start"] and position <= gene["stop"])
	
	def convert_left_locations_to_genes(self,dataset):
		print "Error wrong call; outdated software version!"
		
		import sys
		sys.exit()
		
		"""
		Gene sorting:
		[  ]
		  [    ]
		   [ ]    <- pitfall with loops
		    [ ]   <- pitfall with loops
		    [  ]
		    |
		
		
		1st:
		 scan all genes with GENE-START <= fusion.pos
		 remember pos_1
		2nd:
		 scan UNTIL gene with GENE-START > fusion.pos  (while GENE-START <= fusion.pos
		
		next iteration start from pos_1
		
		
		
		
		
		
		Possibility:
		
		ref:
		      [     ]     <- 1st gene
		 |                <- 1st event
		
		
		ref:
		      [     ]     <- 1st gene
		        |         <- 1st event
		
		ref:
		      [     ]     <- 1st gene
		               |  <- 1st event
		
		"""
		
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
	
	def show_me(self):
		for chromosome in self.annotations_left_indexed.keys():
			print "Chromosome: "+str(chromosome)
			for gene in self.annotations_left_indexed[chromosome]:
				print gene
		print "---------------------"
