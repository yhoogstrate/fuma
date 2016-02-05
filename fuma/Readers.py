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

import logging
import re

from Fusion import Fusion
from FusionDetectionExperiment import FusionDetectionExperiment

class ReadCGhighConfidenceJunctionsBeta(FusionDetectionExperiment):
	logger = logging.getLogger("FuMA::Readers::ReadCGhighConfidenceJunctionsBeta")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		
		self.parse_left_chr_column = -1
		self.parse_right_chr_column = -1
		
		self.parse_left_pos_column = -1
		self.parse_right_pos_column = -1
		
		self.parse_sequence_column = -1
		self.parse_transition_sequence_column = -1
		
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		
		if(len(line) > 0):
			if(line[0] != "#"):
				if(line[0] == ">"):
					self.parse_line__header(line)
				else:
					self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line[1:]
		line = line.split("\t")
		
		self.parse_left_chr_column = line.index("LeftChr")
		self.parse_right_chr_column = line.index("RightChr")
		
		self.parse_left_pos_column = line.index("LeftPosition")
		self.parse_right_pos_column = line.index("RightPosition")
		
		self.parse_sequence_column = line.index("AssembledSequence")
		self.parse_transition_sequence_column = line.index("TransitionSequence")
		
		self.parse_left_strand = line.index("LeftStrand")
		self.parse_right_strand = line.index("RightStrand")
		
		self.parse_id = line.index("Id")
	
	def parse_line__fusion(self,line):
		line = line.split("\t")
		
		left_chr = line[self.parse_left_chr_column]
		right_chr = line[self.parse_right_chr_column]
		left_pos = line[self.parse_left_pos_column]
		right_pos = line[self.parse_right_pos_column]
		
		if(self.parse_sequence_column >= len(line)):
			sequence = ""
		else:
			sequence = line[self.parse_sequence_column]
		
		transition_sequence = line[self.parse_transition_sequence_column]
		left_strand = line[self.parse_left_strand]
		right_strand = line[self.parse_right_strand]
		
		if not left_strand or not right_strand:
			left_strand = None
			right_strand = None
		
		f = Fusion(left_chr, right_chr, left_pos, right_pos, sequence, transition_sequence, left_strand, right_strand,self.name)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadIlluminaHiSeqVCF(FusionDetectionExperiment):
	logger = logging.getLogger("FuMA::Readers::ReadCGhighConfidenceJunctionsBeta")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.breaks = {}
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		if(len(line) > 0):
			if(line[0] != "#"):
				line = line.split("\t")
				
				sv_type = line[7].split("SVTYPE=",1)[1].split(";",1)[0]
				if(sv_type == "DEL"):
					end = line[7].split("END=",1)[1].split(";",1)[0]
					
					f = Fusion(line[0],line[0],line[1],end,False,line[3],None,None,self.name)
					f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[2],'dataset':f.dataset_name})# Secondary location(s)
					
					self.add_fusion(f)
					
				elif(sv_type == "BND"):
					mate = line[7].split("MATEID=",1)[1].split(";",1)[0]
					self.breaks[line[2]] = {'line':line,'mate':mate}
	
	def process_mates(self):
		while(len(self.breaks) > 0):
			item_1 = self.breaks.keys()[0]
			item_2 = self.breaks[item_1]["mate"]
			
			if(self.breaks.has_key(item_2)):
				line_1 = self.breaks[item_1]["line"]
				line_2 = self.breaks[item_2]["line"]
				
				f = Fusion(line_1[0],line_2[0],line_1[1],line_2[1],False,line_1[3],None,None,self.name)
				f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line_1[2],'dataset':f.dataset_name})# Secondary location(s)
				
				self.add_fusion(f)
				
				del(self.breaks[item_2])
				
				f.show_me()
			else:
				self.logger.error("Inappropriate file - missing link to: "+item_2)
			
			del(self.breaks[item_1])
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.i = 0
		with open(self.filename,"r") as fh:
			for line in fh:
				self.i += 1
				self.parse_line(line)
		
		self.process_mates()
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadTophatFusionPre(FusionDetectionExperiment):
	"""Parses Tophat Fusion's file 'fusions.out'
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadTophatFusionPre")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		if(len(line) > 0):
			q = line
			line = line.split("@")
			line[0] = line[0].split("\t")
			line[2] = line[2].strip().split(" ")
			line[3] = line[3].strip().split(" ")
			
			chromosomes = line[0][0].split("-")
			
			if(len(line[3]) > 1):
				sequence = line[2][0] + line[3][1]
			else:
				sequence = False
			
			f = Fusion(chromosomes[0],chromosomes[1],line[0][1],line[0][2],sequence,False,line[0][3][0],line[0][3][1],self.name)
			f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
			
			self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.i = 0
		with open(self.filename,"r") as fh:
			for line in fh:
				self.i += 1
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadTophatFusionPostPotentialFusion(FusionDetectionExperiment):
	"""Parsess TopHat Fusion post's output file 'potential_fusion.txt'
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadTophatFusionPostPotentialFusion")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		
		self.parse()
	
	def reset(self):
		self.chr_1 = False
		self.chr_2 = False
		
		self.break_1 = False
		self.break_2 = False
		
		self.seq = False
		self.insert_seq = False
	
	def parse_line_type_0(self,line):
		#self._id = line.strip()										# this is the sample-id; tophat expects multiple experiments per tophat-fusion post run
		
		line = line.strip().split(" ")
		_chr = line[1].split("-")
		
		self.chr_1 = _chr[0]
		self.chr_2 = _chr[1]
		
		self.break_1 = int(line[2])
		self.break_2 = int(line[3])
		
		self.left_strand = line[4][0]
		self.right_strand = line[4][1]
	
	def parse_line_type_1(self,line):
		line = line.strip().split(" ")
		if(len(line) > 0):
			self.seq = line[0]
		
	def parse_line_type_2(self,line):
		line = line.strip().split(" ")
		if(len(line) > 1):
			self.seq += line[1]
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.reset()
		
		with open(self.filename,"r") as fh:
			i = 0
			j = 0
			for line in fh:
				line_type = i % 6
				if(line_type == 0):
					self.parse_line_type_0(line)
					j += 1
				elif(line_type == 1):
					self.parse_line_type_1(line)
				elif(line_type == 2):
					self.parse_line_type_2(line)
					
					f = Fusion(self.chr_1,self.chr_2,self.break_1,self.break_2,self.seq,self.insert_seq,self.left_strand,self.right_strand,self.name)
					f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(j),'dataset':f.dataset_name})# Secondary location(s)
					
					self.add_fusion(f)
				
				i += 1
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadTophatFusionPostResult(FusionDetectionExperiment):
	"""Parsess TopHat Fusion post's output file 'result.txt'
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadTophatFusionPostResult")
	
	parse_left_chr_column = 2
	parse_right_chr_column = 5
	
	break_left = 3
	break_right = 6
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.i = 0
		
		self.filename = arg_filename
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line_stripped = line.strip()
		if(len(line) > 0):
			line = line.split("\t")
			
			f = Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],line[self.break_left],line[self.break_right],None,False,None,None,self.name)
			f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})
			
			self.i += 1
			
			self.add_fusion(f)



class ReadTophatFusionPostResultHtml(FusionDetectionExperiment):
	"""Parsess TopHat Fusion post's output file 'result.html'
	
	Regex lines similar to:
<P><P><P><BR>
2. chr22-chr9 ff
<TABLE CELLPADDING=3 BORDER="1">
<TR><TD ALIGN="LEFT"><a href="#fusion_7">sample_1</a></TD>
<TD ALIGN="LEFT">NASP</TD>
<TD ALIGN="LEFT">chr1</TD>
<TD ALIGN="RIGHT">46070687</TD>
<TD ALIGN="LEFT">ENSG00000254777</TD>
<TD ALIGN="LEFT">chr8</TD>
<TD ALIGN="RIGHT">61852322</TD>
<TD ALIGN="RIGHT"><a href="#read_7">11</a></TD>
<TD ALIGN="RIGHT"><a href="#pair_7">424</a></TD>
<TD ALIGN="RIGHT">1</TD>
</TR>
</TABLE>
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadTophatFusionPostResultHtml")
	
	parse_left_chr_column = 2
	parse_right_chr_column = 5
	
	break_left = 3
	break_right = 6
	
	ppp_block_match = re.compile('<P><P><P><BR>[^0-9]+[0-9]+\. [^\n]*?(.)(.)\n(<TABLE .*?</TABLE>)',re.S)
	
	td_match = ".*?<TD [^>]+>([^<]+)</TD>"
	table_block_match = re.compile('href="#fusion_([^"]+)">'+td_match+td_match+td_match+td_match+td_match+td_match,re.S)
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		fh = open(self.filename,"r")
		content = fh.read()
		fh.close()
		
		self.parse_ppp_blocks(content)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_ppp_blocks(self,content):
		for match in re.finditer(self.ppp_block_match,content):
			match = match.groups()
			self.parse_table_blocks(match[0],match[1],match[2])
	
	def parse_table_blocks(self,strand1,strand2,ppp_block):
		for match in re.finditer(self.table_block_match,ppp_block):
			match = match.groups()
			
			f = Fusion(match[2] ,match[5] ,match[3] ,match[6] ,None,False,strand1,strand2,self.name)
			f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':match[0],'dataset':f.dataset_name})
			
			self.add_fusion(f)



class ReadDefuse(FusionDetectionExperiment):
	"""
	Defuse's gene annotations are DEFINITELY 0-based (both start as end)
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadDefuse")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		if(self.parse_header):
			self.parse_line__header(line)
		else:
			self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line.strip().split("\t")
		
		self.parse_left_chr_column = line.index("gene_chromosome1")
		self.parse_right_chr_column = line.index("gene_chromosome2")
		
		self.parse_left_pos_column = line.index("genomic_break_pos1")
		self.parse_right_pos_column = line.index("genomic_break_pos2")
		
		self.parse_sequence_column = line.index("splitr_sequence")
		
		self.parse_left_strand_column = line.index("genomic_strand1")
		self.parse_right_strand_column = line.index("genomic_strand2")
		
		self.parse_id = line.index("cluster_id")
		
		self.parse_header = False
	
	def parse_line__fusion(self,line):
		line = line.strip().split("\t")
		
		left_pos = int(line[self.parse_left_pos_column])
		right_pos = int(line[self.parse_right_pos_column])
		
		f = Fusion( \
			line[self.parse_left_chr_column], \
			line[self.parse_right_chr_column], \
			left_pos, \
			right_pos, \
			line[self.parse_sequence_column], \
			False, \
			line[self.parse_left_strand_column], \
			line[self.parse_right_strand_column], \
			self.name \
		)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadFusionMap(FusionDetectionExperiment):
	logger = logging.getLogger("FuMA::Readers::ReadFusionMap")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse_header = False
		
		self.parse()
	
	def parse_line__header(self,line):
		line = line.strip().split("\t")
		
		self.parse_left_chr_column = line.index("Chromosome1")
		self.parse_right_chr_column = line.index("Chromosome2")
		
		self.parse_left_pos_column = line.index("Position1")
		self.parse_right_pos_column = line.index("Position2")
		
		self.parse_strand_columns = line.index("Strand")
		
		self.parse_id_column = line.index("FusionID")
		
		self.parse_header = True
	
	def parse_line__fusion(self,line):
		line = line.strip().split("\t")
		
		f = Fusion(
			line[self.parse_left_chr_column], \
			line[self.parse_right_chr_column], \
			line[self.parse_left_pos_column], \
			line[self.parse_right_pos_column], \
			False, \
			False, \
			line[self.parse_strand_columns][0], \
			line[self.parse_strand_columns][1], \
			self.name \
		)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[self.parse_id_column],'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					if(self.parse_header == False):
						self.parse_line__header(line)
					else:
						self.parse_line__fusion(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadChimeraScanAbsoluteBEDPE(FusionDetectionExperiment):
	"""
		ChimeraScan provides two types of BEDPE files. We classify them
		in the "absolute" and the "relative" BEDPE files:
		- The "absolute" files have absolute genomic coordinates, relative
		to the chromosome.
		- The "relative" files contrain coordinates relative to genes,
		unless they are not in gene regions of course. To convert these
		files back to absolute you MUST have the correct gene annotation
		file. Therefore this parser is not suited, since it only needs
		one file as input.
		
		Chimerascan provides no clear definition on how to extract
		exact breakpoints from the output files. The BEDPE files contain
		the columns 'start5p', 'end5p', 'start3p' and 'end3p' where we
		observe that always ('start5p' < 'end5p') and
		('start3p' < 'end3p'). Therefore we assume they are the regions
		that are covered by reads. By using use the columns 'strand5p'
		and 'strand3' we hope to find the exact breakpoints:
		
		if 'strand5p' == '+'
			breakpoint_1 = end5p
		elif 'strand5p' == '-'
			breakpoint_1 = start5p
		
		if 'strand3p' == '+'
			breakpoint_2 = end3p
		elif 'strand3p' == '-'
			breakpoint_2 = start3p
		
		Probably for those reads that have a non-empty column
		'breakpoint_spanning_reads' it is possible to estimate the
		breakpoint more in-depth, but we don't know how to implement
		this (yet?).
		
		@todo ensure that this is also applied in the conversion file(s)
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadChimeraScanAbsoluteBEDPE")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		if(self.parse_header):
			if(line[0] == "#"):
				self.parse_line__header(line)
			else:
				self.parse_left_chr_column = 0
				self.parse_right_chr_column = 3
				
				self.parse_start5p_column = 1
				self.parse_end5p_column = 2
				
				self.parse_start3p_column = 4
				self.parse_end3p_column = 5
				
				self.parse_left_strand = 8
				self.parse_right_strand = 9
				
				self.parse_id = 6
				
				self.parse_line__fusion(line)
			
			self.parse_header = False
		else:
			self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line.strip().split("\t")
		
		self.parse_left_chr_column = line.index("#chrom5p")
		self.parse_right_chr_column = line.index("chrom3p")
		
		self.parse_start5p_column = line.index("start5p")
		self.parse_end5p_column = line.index("end5p")
		
		self.parse_start3p_column = line.index("start3p")
		self.parse_end3p_column = line.index("end3p")
		
		self.parse_left_strand = line.index("strand5p")
		self.parse_right_strand = line.index("strand3p")
		
		self.parse_id = line.index("chimera_cluster_id")
	
	def parse_line__fusion(self,line):
		line = line.strip().split("\t")
		
		if(line[self.parse_left_strand] == "+"):
			breakpoint_1 = int(line[self.parse_end5p_column])-1			# BEDPE end-positions are 1-based: http://bedtools.readthedocs.org/en/latest/content/general-usage.html#bedpe-format
		elif(line[self.parse_left_strand] == "-"):
			breakpoint_1 = line[self.parse_start5p_column]
		
		if(line[self.parse_right_strand] == "+"):
			breakpoint_2 = int(line[self.parse_end3p_column])-1			# BEDPE end-positions are 1-based: http://bedtools.readthedocs.org/en/latest/content/general-usage.html#bedpe-format
		elif(line[self.parse_right_strand] == "-"):
			breakpoint_2 = line[self.parse_start3p_column]
		
		f = Fusion( \
			line[self.parse_left_chr_column], \
			line[self.parse_right_chr_column], \
			breakpoint_1, \
			breakpoint_2, \
			False, \
			False, \
			line[self.parse_left_strand], \
			line[self.parse_right_strand], \
			self.name \
		)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def convert_to_absolute_coordinates(self,gene_features,output):
		for fusion in self:
			gene_id_left = fusion.get_left_chromosome()#[3:]
			gene_id_right = fusion.get_right_chromosome()#[3:]
			
			go=[]
			try:
				left_gene = gene_features.index[gene_id_left]
			except:
				go.append(gene_id_left)
			
			try:
				right_gene = gene_features.index[gene_id_right]
			except:
				go.append(gene_id_right)
			
			if(len(go) == 0):
				new_left_pos = left_gene[1] + fusion.get_left_break_position() - 1
				new_right_pos = right_gene[1] + fusion.get_right_break_position() - 1
				
				fusion.set(left_gene[0],right_gene[0],new_left_pos,new_right_pos,fusion.sequence,fusion.transition_sequence,fusion.left_strand,fusion.right_strand)
			else:
				self.logger.warning("Can not find genes with id: "+", ".join(go))
		
		self.export_to_CG_Junctions_file(output)



class ReadFusionCatcherFinalList(FusionDetectionExperiment):
	"""Reads the FusionCatchers 'final-list_candidate-fusion-genes.txt'
	"""
	
	logger = logging.getLogger("FuMA::Readers::ReadFusionCatcherFinalList")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		
		if(len(line) > 0):
			if(self.parse_header):
				self.parse_line__header(line)
			else:
				self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line.split("\t")
		
		try:
			self.parse_left_column = line.index("Fusion_point_for_gene_1(5end_fusion_partner)")
			self.parse_right_column = line.index("Fusion_point_for_gene_2(3end_fusion_partner)")
		except:
			self.parse_left_column = line.index("Fusion_gene_1_position(5end_partner)")
			self.parse_right_column = line.index("Fusion_gene_2_position(3end_partner)")
		
		self.parse_sequence_column = line.index("Fusion_sequence")
		
		self.parse_header = False
	
	def parse_line__fusion(self,line):
		line = line.split("\t")
		
		left_chr,left_pos,left_strand = line[self.parse_left_column].split(":")
		right_chr,right_pos,right_strand = line[self.parse_right_column].split(":")
		
		f = Fusion(
			left_chr, \
			right_chr, \
			int(left_pos), \
			int(right_pos), \
			line[self.parse_sequence_column], \
			False, \
			left_strand, \
			right_strand, \
			self.name \
			)
		
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
		
		self.i += 1
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		self.i = 1
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))


class FusionCatcherIndices:
	logger = logging.getLogger("FuMA::Readers::FusionCatcherIndices")
	
	def __init__(self):
		self.gene_index = {}
		self.transcript_index = {}
		self.exon_index = {}
	
	def parse_genes(self,filename):
		self.logger.info("Parsing file (for genes): "+str(filename))
		
		i = 0
		with open(filename,"r") as fh:
			for line in fh:
				line_s = line.strip("\n")
				if(len(line_s) > 0):
					params = line_s.split("\t")
					
					gene_id = params[0]
					start = int(params[1])
					stop = int(params[2])
					
					self.gene_index[gene_id] = {'chromosome': "chr"+params[4],'start': min(start,stop), 'stop': max(start,stop)}
					self.gene_index[gene_id]['center'] = int(round(0.5*(self.gene_index[params[0]]['start'] + self.gene_index[params[0]]['stop'] + 1)))
					
					i += 1
		
		self.logger.info("Parsed genes: "+str(i))
	
	def parse_transcripts_line(self,line):
		start = None
		stop = -1
		chromosome = 'unknown'
		for entry in line.split(";"):
			if(entry[0:3] == "ex="):
				for param in entry.split(","):
					key,value = param.split("=")
					if(key == "sc"):
						value = int(value)
						if(start == None or value < start):
							start = value
					elif(key == "ec"):
						value = int(value)
						if(value > stop):
							stop = value
			elif(entry[0:4] == "chr="):
				chromosome = "chr"+entry.split("=")[1]
		return {'chromosome': chromosome, 'start': start,'stop': stop}
	
	def parse_transcripts(self,filename):
		self.logger.info("Parsing file (for transcripts): "+str(filename))
		
		i = 0
		
		with open(filename,"r") as fh:
			for line in fh:
				line_s = line.strip("\n")
				if(len(line_s) > 0):
					params = line_s.split("\t")
					transcript_id = params[0].split(";")[0]
					self.transcript_index[transcript_id] = self.parse_transcripts_line(line_s)
					self.transcript_index[transcript_id]['center'] = int(round(0.5*(self.transcript_index[transcript_id]['start'] + self.transcript_index[transcript_id]['stop'] + 1)))
					
					i += 1
		
		self.logger.info("Parsed transcripts: "+str(i))
	
	def parse_exons(self,filename):
		"""Exons file is of following syntax:
		Pathway-id \t Gene-id \t Transcript-id \t Exon-id \t start \t stop .... \t chromosome
		"""
		self.logger.info("Parsing file (for exons): "+str(filename))
		
		with open(filename,"r") as fh:
			for line in fh:
				line_s = line.strip("\n")
				
				if(len(line_s) > 0):
					params = line_s.split("\t")
					
					exon_id = params[3]
					start = int(params[4])
					stop = int(params[5])
					
					self.exon_index[exon_id] = {'chromosome': "chr"+params[12],'start': min(start,stop), 'stop': max(start,stop)}
					self.exon_index[exon_id]['center'] = int(round(0.5*(self.exon_index[exon_id]['start'] + self.exon_index[exon_id]['stop'] + 1)))



class ReadFusionCatcherMAP(FusionDetectionExperiment):
	logger = logging.getLogger("FuMA::Readers::ReadFusionCatcherMAP")
	
	def __init__(self,arg_filename,name,references):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.references = references
		
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		self.i = 1
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line = line.strip("\n")
		if(len(line) > 0):
			params = line.split("\t")
			items = params[2].split(";")
			exons = items[2].split("-")
			
			exon1 = self.references.exon_index[exons[0]]
			exon2 = self.references.exon_index[exons[1]]
			
			f = Fusion(exon1['chromosome'],exon2['chromosome'],exon1['center'],exon2['center'],None,False,None,None,self.name)
			f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
			
			self.add_fusion(f)
			
			self.i += 1



class ReadFusionCatcherPreliminaryList(FusionDetectionExperiment):
	logger = logging.getLogger("FuMA::Readers::ReadFusionCatcherPreliminaryList")
	
	parse_left_gene = 0
	parse_right_gene = 1
	
	def __init__(self,arg_filename,name,references):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.references = references
		
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.parse_header = True
		self.i = 0
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line = line.strip("\n")
		if(len(line) > 0):
			if(self.i >= 1):
				params = line.split("\t")
				
				gene1 = self.references.gene_index[params[self.parse_left_gene]]
				gene2 = self.references.gene_index[params[self.parse_right_gene]]
				
				f = Fusion(gene1['chromosome'],gene2['chromosome'],gene1['center'],gene2['center'],None,False,None,None,self.name)
				f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
				
				self.add_fusion(f)
			
			self.i += 1


class ReadRNASTARChimeric(FusionDetectionExperiment):
	"""Example file syntax:
chr8	70572329	+	chr8	70572307	-	0	0	1	HWI-1KL113:71:D1G2NACXX:1:1102:12932:160607	70572289	40M61S	70572146	101M-1p61M40S
chr8	29921084	-	chr8	29921059	+	0	0	0	HWI-1KL113:71:D1G2NACXX:1:1102:16721:160648	29921085	67S34M-11p101M	29921060	34S67M
chr7	99638140	+	chr7	99638098	-	0	0	3	HWI-1KL113:71:D1G2NACXX:1:1102:17025:160706	99637628	44M1I56M364p48M53S	99638045	53M48S
	"""
	parse_left_chr_column = 0
	parse_left_pos_column = 1
	parse_left_strand_column = 2
	
	parse_right_chr_column = 3
	parse_right_pos_column = 4
	parse_right_strand_column = 5
	
	logger = logging.getLogger("FuMA::Readers::ReadRNASTARChimeric")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.i = 1
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
					self.i += 1
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line = line.strip().split("\t")
		
		left_pos = int(line[self.parse_left_pos_column])
		right_pos = int(line[self.parse_right_pos_column])
		
		f = Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],left_pos,right_pos,None,False,line[self.parse_left_strand_column],line[self.parse_right_strand_column],self.name)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)


class ReadRNASTARFusionFinal(FusionDetectionExperiment):
	"""Example file syntax:
#fusion_name	JunctionReads	SpanningFrags	Splice_type	LeftGene	LeftBreakpoint	RightGene	RightBreakpoint	JunctionReads	SpanningFrags
SLC12A7--AAMDC	13	55	ONLY_REF_SPLICE	SLC12A7^ENSG00000113504.15	chr5:1085347:-	AAMDC^ENSG00000087884.10	chr11:77580768:+	SRR018266.9510307,SRR018266.322525,SRR018266.12147603,SRR018266.7829232,SRR018266.4364047,SRR018266.1580199,SRR018266.5313966,SRR018266.12643509,SRR018266.3302049,SRR018266.940865,SRR018266.3535489,SRR018266.1693248,SRR018266.3380501	SRR018266.7168134,SRR018266.10866101,SRR018266.3328842,SRR018266.4816578,SRR018266.13389863,SRR018266.9676001,SRR018266.8967505,SRR018266.5967402,SRR018266.13331054,SRR018266.11445954,SRR018266.2955213,SRR018266.1926608,SRR018266.1756809,SRR018266.9216697,SRR018266.664635,SRR018266.1694889,SRR018266.13966627,SRR018266.702434,SRR018266.14400502,SRR018266.9945013,SRR018266.7196489,SRR018266.9715412,SRR018266.5864824,SRR018266.4363568,SRR018266.846284,SRR018266.5640930,SRR018266.8618635,SRR018266.3614974,SRR018266.1146973,SRR018266.5240308,SRR018266.3337293,SRR018266.12915124,SRR018266.10757122,SRR018266.407787,SRR018266.11176512,SRR018266.6878929,SRR018266.3867682,SRR018266.6620348,SRR018266.14553851,SRR018266.11006496,SRR018266.4887599,SRR018266.4374589,SRR018266.13353188,SRR018266.12529837,SRR018266.8140358,SRR018266.9425841,SRR018266.3425026,SRR018266.434919,SRR018266.10983920,SRR018266.8446368,SRR018266.7921877,SRR018266.12225667,SRR018266.8457671,SRR018266.1613806,SRR018266.4055216	
CCT3--C1orf61	6	42	ONLY_REF_SPLICE	CCT3^ENSG00000163468.10	chr1:156294763:-	C1orf61^ENSG00000125462.12	chr1:156374393:-	SRR018266.9049684,SRR018266.13478646,SRR018266.3340105,SRR018266.14621241,SRR018266.11203077,SRR018266.6789091	SRR018266.9000529,SRR018266.396588,SRR018266.2731797,SRR018266.10679985,SRR018266.8373243,SRR018266.2699994,SRR018266.10010923,SRR018266.385212,SRR018266.11975157,SRR018266.7529212,SRR018266.6838575,SRR018266.14665000,SRR018266.10947416,SRR018266.10887397,SRR018266.13919642,SRR018266.7053633,SRR018266.10738050,SRR018266.6926730,SRR018266.4054306,SRR018266.8555211,SRR018266.10999663,SRR018266.11295642,SRR018266.11806376,SRR018266.11254990,SRR018266.6631020,SRR018266.5538318,SRR018266.2792602,SRR018266.7341944,SRR018266.8837183,SRR018266.10661720,SRR018266.6529302,SRR018266.1168740,SRR018266.7593177,SRR018266.7342096,SRR018266.8966067,SRR018266.10860075,SRR018266.12337779,SRR018266.4494592,SRR018266.14709847,SRR018266.7909211,SRR018266.3892227,SRR018266.11524233	
EPB41--YIPF3	5	0	INCL_NON_REF_SPLICE	EPB41^ENSG00000159023.14	chr1:29446010:+	YIPF3^ENSG00000137207.7	chr6:43479658:-	SRR018266.13578666,SRR018266.12258196,SRR018266.5607190,SRR018266.14197979,SRR018266.6975350	.	
	"""
	logger = logging.getLogger("FuMA::Readers::ReadRNASTARChimeric")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.header = None
		
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		self.i = 1
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					if(self.header):
						self.parse_line(line)
						self.i += 1
					else:
						self.parse_header_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line = line.strip().split("\t")
		
		left_break = line[self.header["parse_left_column"]].split(":")
		right_break = line[self.header["parse_right_column"]].split(":")
		
		f = Fusion(
					left_break[0], \
					right_break[0], \
					left_break[1], \
					right_break[1], \
					None, \
					False, \
					left_break[2], \
					right_break[2], \
					self.name)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':str(self.i),'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse_header_line(self,line):
		line = line.split("\t")
		
		self.header = {}
		self.header["parse_left_column"] = line.index("LeftBreakpoint")
		self.header["parse_right_column"] = line.index("RightBreakpoint")



class ReadChimeraPrettyPrint(FusionDetectionExperiment):
	""" Example file sytax:
"gene1"	"chr.gene1"	"breakpoint.gene1"	"strand.gene1"	"transcripts.gene1"	"gene2"	"chr.gene2"	"breakpoint.gene2"	"strand.gene2"	"transcripts.gene2"	"fusion.breakpoint"	"supporting.reads"
"MT-ND5"	"chrMT"	"14006"	"+"	"NA"	"J01415.25"	"chrMT"	"8407"	"-"	"NA"	"TAAAATAAAATCCCC"	"11"
"MT-ND4"	"chrMT"	"11706"	"-"	"NA"	"MT-ND2"	"chrMT"	"5320"	"+"	"NA"	"GTATAATACGCCTTC"	"99"
	"""
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		self.filename = arg_filename
		self.columns = None
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line.strip())
	
	def parse_line(self,line):
		if len(line) > 0:
			if self.columns == None:
				self.parse_line__header(line)
			else:
				self.parse_line__fusion(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def cleanup_params(self,params):
		params_clean = []
		
		for param in params:
			if len(param) > 1 and param[0] == '"' and param[-1] == '"':
				params_clean.append(param[1:-1].replace('\\"','"'))
			else:
				params_clean.append(param)
		
		return params_clean
	
	def parse_line__header(self,line):
		params = self.cleanup_params(line.split("\t"))
		
		self.columns = { \
			'left_chr':     params.index('chr.gene1') , \
			'right_chr':    params.index('chr.gene2') , \
			'left_pos':     params.index('breakpoint.gene1') , \
			'right_pos':    params.index('breakpoint.gene2') , \
			'left_strand':  params.index('strand.gene1') , \
			'right_strand': params.index('strand.gene2') , \
			'transequence': params.index('fusion.breakpoint') }
	
	def parse_line__fusion(self,line):
		line = self.cleanup_params(line.split("\t"))
		
		left_chr = line[self.columns['left_chr']]
		right_chr = line[self.columns['right_chr']]
		
		left_pos = line[self.columns['left_pos']]
		right_pos = line[self.columns['right_pos']]
		
		left_strand = line[self.columns['left_strand']]
		right_strand = line[self.columns['right_strand']]
		
		if(self.columns['transequence'] != "NA"):
			transition_sequence = line[self.columns['transequence']]
		else:
			transition_sequence = None
		
		uid = str(len(self))
		
		f = Fusion(left_chr, right_chr, left_pos, right_pos, None, transition_sequence, left_strand, right_strand, self.name)
		f.add_location({ \
			'left':[f.get_left_chromosome(), \
			f.get_left_break_position()], \
			'right':[f.get_right_chromosome(), \
			f.get_right_break_position()], \
			'id':uid, \
			'dataset':f.dataset_name \
		})# Secondary location(s)
		
		self.add_fusion(f)



class Read123SVDeNovo(FusionDetectionExperiment):
	"""
	Example file syntax:
1	157121267	157128069	1	157778480	157785127	6(108)	TT(46)tt(62)	108	100	-122	0.177571632724291	inversion
1	157778113	157784838	1	178271911	178279982	6(64)	HH(32)hh(32)	64	100	-478	0.156719992865115	inversion
9	114014626	114014651	9	114014626	114014654	6(59)	TH(37)th(22)	59	100	-26	0.42919821329749	insertion
	"""
	
	parse_left_chr_column = 0
	parse_right_chr_column = 3
	
	parse_left_pos_column = [1,2]
	parse_right_pos_column = [4,5]
	
	logger = logging.getLogger("FuMA::Readers::Read123SVDeNovo")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		
		self.parse()
	
	def parse_line(self,line):
		line = line.strip()
		
		if(len(line) > 0):
			if(line[0] != "#"):
				self.parse_line__fusion(line)
	
	def parse_line__fusion(self,line):
		line = line.split("\t")
		
		left_chr = line[self.parse_left_chr_column]
		right_chr = line[self.parse_right_chr_column]
		
		"""
		Columns 1-3: Chromosome, start position and end position of range where first breakpoint can be found.
		Columns 4-6: Chromosome, start position and end position of range where second breakpoint can be found.
		
		for left and right position, pick the mean of both...
		"""
		
		left_pos = (int(line[self.parse_left_pos_column[0]]) + int(line[self.parse_left_pos_column[1]])) / 2
		right_pos = (int(line[self.parse_right_pos_column[0]]) + int(line[self.parse_right_pos_column[1]])) / 2
		
		sequence = ""
		
		transition_sequence = ""
		left_strand = "+"
		right_strand = "+"
		
		"""
		@todo: Column 8: Relative orientation and strand of the tag pairs in the SV call. Specifies the number of observed tags that link first segment to second in "tail-to-head" (TH, 3-prime of first segment is linked to 5-prime of second segment) or any other combination. Strand of first read in the pair is reflected as well: big caps for plus strand (e.g. 'TH') and small caps for minus strand (e.g. 'th').
		@link http://tools.genomes.nl/123sv.html
		"""
		
		uid = line[0]+":"+line[1]+","+line[2]+"-"+line[3]+":"+line[4]+","+line[5]
		
		f = Fusion(left_chr, right_chr, left_pos, right_pos, sequence, transition_sequence, left_strand, right_strand,self.name)
		f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':uid,'dataset':f.dataset_name})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))





class ReadOncofuse(FusionDetectionExperiment):
	"""
	Example file syntax:
SAMPLE_ID	FUSION_ID	TISSUE	SPANNING_READS	ENCOMPASSING_READS	GENOMIC	5_FPG_GENE_NAME	5_IN_CDS?	5_SEGMENT_TYPE	5_SEGMENT_ID	5_COORD_IN_SEGMENT	5_FULL_AA	5_FRAME	3_FPG_GENE_NAME	3_IN_CDS?	3_SEGMENT_TYPE	3_SEGMENT_ID	3_COORD_IN_SEGMENT	3_FULL_AA	3_FRAME	FPG_FRAME_DIFFERENCE	P_VAL_CORR	DRIVER_PROB	EXPRESSION_GAIN	5_DOMAINS_RETAINED	3_DOMAINS_RETAINED	5_DOMAINS_BROKEN	3_DOMAINS_BROKEN	5_PII_RETAINED	3_PII_RETAINED	CTF	G	H	K	P	TF
Chimeric.out.junction	124	EPI	3	0	chr19:58421250>chr19:58370230	ZNF417	Yes	Exon	3	232	131	1	ZNF587	Yes	Exon	3	287	427	1	2	0.17279408624033415	0.9794409302490255	-0.5987947056878361	Krueppel-associated box[Domain];Zinc finger, C2H2-like[Domain];Zinc finger, C2H2[Domain]	Zinc finger C2H2-type/integrase DNA-binding domain[Domain];Zinc finger, C2H2-like[Domain];Zinc finger, C2H2[Domain]	Zinc finger, C2H2[Domain]				0.02895169933768111	0.0	7.449265196664073E-5	0.0	0.005733585138730341	0.09019563359624447
Chimeric.out.junction	136	EPI	2	0	chr2:27293455>chr2:27293539	AGBL5	No	Exon	15	496	887	0	OST4	No	Exon	3	50	0	0	0	1.0	6.333147056899562E-4	NaN	Peptidase M14, carboxypeptidase A[Domain]						0.0	0.0	0.0	0.0	0.0	0.0
Chimeric.out.junction	183	EPI	2	0	chr3:52300997>chr19:36726707	WDR82	Yes	Exon	3	37	98	1	ZNF146	No	Exon	3	147	293	0	1	0.09094039431674192	0.9933018406452259	0.05781956186247239	WD40 repeat[Repeat];G-protein beta WD-40 repeat[Repeat]	Zinc finger C2H2-type/integrase DNA-binding domain[Domain];Zinc finger, C2H2-like[Domain];Zinc finger, C2H2[Domain]	WD40-repeat-containing domain[Domain];WD40 repeat[Repeat];WD40/YVTN repeat-like-containing domain[Domain]				0.0159434572349889	0.0	3.7249030480169165E-5	0.0	0.0034150577334664174	0.05079275660405341
	"""
	
	parse_genomic_column = 5
	parse_fusionid_column = 1
	
	logger = logging.getLogger("FuMA::Readers::ReadOncofuse")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		self.filename = arg_filename
		self.parse()
	
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line = line.strip()
				if(len(line) > 0):
					self.parse_line(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_line(self,line):
		line = line.strip().split("\t")
		
		if(line[0] != "SAMPLE_ID"):										# Skip header line
			genomic = line[self.parse_genomic_column]
			left,right = genomic.split(">")
			left = left.split(":")
			right = right.split(":")
			
			f = Fusion(left[0],right[0],int(left[1]),int(right[1]),None,False,None,None,self.name)
			f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':line[self.parse_fusionid_column],'dataset':f.dataset_name})# Secondary location(s)
			
			self.add_fusion(f)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))



class ReadTrinityGMAP(FusionDetectionExperiment):
	regexes = {}
	#regexes['Path'] = 'Path [12]: query ([\S]*?)\.\.([\S]*?) \(([0-9\-]+) bp\) [\S]+ [\S]+ ([^: ]+):([^\.]+)\.\.([^ ]+) \(([0-9\-]+) bp\)'
	#regexes['cDNA direction'] = False
	regexes['Genomic pos'] = 'Genomic pos: ([^:]+):([^\.]+)\.\.([^ ]+) \(([^ ]+) strand\)'
	regexes['Accessions'] = 'Accessions: ([^:]+):([^\.]+)\.\.([^ ]+) \(out of ([^ ]+) bp\)'
	#regexes['Number of exons'] = 'Number of exons: ([0-9]+)'
	#regexes['Coverage'] = 'Coverage: ([^ ]+) \(query length: ([^ ]+) bp\)'
	#regexes['Trimmed coverage'] = 'Trimmed coverage: ([^ ]+) \(trimmed length: ([^ ]+) bp, trimmed region: ([^ ]+)\.\.([^\)]+)\)'
	#regexes['Percent identity'] = 'Percent identity: ([^ ]+) \(([^ ]+) matches, ([^ ]+) mismatches, ([^ ]+) indels, ([^ ]+) unknowns\)'
	#regexes['Translation'] = 'Translation: ([^\.]+)\.\.([^ ]+) \(([^ ]+) aa\)'
	#regexes['Amino acid changes'] = 'Amino acid changes: ([^\n]*)'
	#regexes['Non-intron gaps'] = 'Non-intron gaps: ([^ ]+) openings, ([^ ]+) bases in cdna; ([^ ]+) openings, ([^ ]+) bases in genome'
	
	logger = logging.getLogger("FuMA::Readers::ReadTrinityGMAP")
	
	def __init__(self,arg_filename,name):
		FusionDetectionExperiment.__init__(self,name)
		
		
		self.filename = arg_filename
		self.parse()
		
	def parse(self):
		self.logger.info("Parsing file: "+str(self.filename))
		
		contig_chunk = []
		contig_name = False
		
		with open(self.filename,"r") as fh:
			for line in fh:
				line_stripped = line.strip()
				if(line_stripped != ""):
					if(line[0:5] == ">comp"):
						contig_chunk = []
						contig_name = line
					elif(line[0:11] == "Alignments:"):
						data = self.parse_contig(contig_name,contig_chunk)
						
						if(data[1]["Genomic pos"][3] == "+"):
							left_pos = data[1]["Accessions"][2]
						else:
							left_pos = data[1]["Accessions"][1]
						
						if(data[2]["Genomic pos"][3] == "+"):
							right_pos = data[2]["Accessions"][1]
						else:
							right_pos = data[2]["Accessions"][2]
						
						uid = contig_name.split(" ")[0].lstrip(">")
						
						f = Fusion(data[1]["Accessions"][0],data[2]["Accessions"][0],left_pos,right_pos,False,False,data[1]["Genomic pos"][3],data[2]["Genomic pos"][3],self.name)
						f.add_location({'left':[f.get_left_chromosome(),f.get_left_break_position()],'right':[f.get_right_chromosome(),f.get_right_break_position()],'id':uid,'dataset':f.dataset_name})
						
						distance = f.get_distance()
						if(distance > 100000 or distance == -1):
							self.add_fusion(f)
					else:
						contig_chunk.append(line)
		
		self.logger.info("Parsed fusion genes: "+str(len(self)))
	
	def parse_contig(self,contig_name,contig_chunk):
		path = 0
		paths = {1:[],2:[]}
		
		i = 0
		
		while i < len(contig_chunk) and path < 3:
			line = contig_chunk[i]
			sline = line.lstrip()
			if(sline[0:7] == 'Path 1:'):
				path = 1
			elif(sline[0:7] == 'Path 2:'):
				path = 2
			elif(sline[0:11] == 'Alignments:'):
				path = 3
			
			if(path > 0 and path < 3):
				paths[path].append(sline)
			
			i += 1
		
		for path in paths.keys():
			paths[path] = self.parse_path(paths[path])
		
		paths['name'] = contig_name.strip()
		
		return paths
	
	def parse_path(self,path_chunk):
		keys = {}
		for line in path_chunk:
			key = line.split(': ')
			key = key[0].replace('Path 1','Path').replace('Path 2','Path')
			if(self.regexes.has_key(key)):
				m = re.search(self.regexes[key],line)
				keys[key] = m.groups()
		return keys
