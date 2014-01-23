#!/usr/bin/env python



import re



from Fusion import Fusion
from HighThroughputFusionDetectionExperiment import HighThroughputFusionDetectionExperiment



class ReadCGhighConfidenceJunctionsBeta(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"DNA")
		
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
		
		f = Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],line[self.parse_left_pos_column],line[self.parse_right_pos_column],line[self.parse_sequence_column],line[self.parse_transition_sequence_column],line[self.parse_left_strand],line[self.parse_right_strand],self.name)
		f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.get_dataset_name()})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)



class ReadIlluminaHiSeqVCF(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
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
					
					f = Fusion(line[0],line[0],line[1],end,False,line[3],"+","+",self.name)
					f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':line[2],'dataset':f.get_dataset_name()})# Secondary location(s)
					
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
				
				f = Fusion(line_1[0],line_2[0],line_1[1],line_2[1],False,line_1[3],"+","+",self.name)
				f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':line_1[2],'dataset':f.get_dataset_name()})# Secondary location(s)
				
				self.add_fusion(f)
				
				del(self.breaks[item_2])
				
				f.show_me()
			else:
				print "ERROR: Inappropriate file - missing link to: "+item_2
			
			del(self.breaks[item_1])
	
	def parse(self):
		self.i = 0
		with open(self.filename,"r") as fh:
			for line in fh:
				self.i += 1
				self.parse_line(line)
		self.process_mates()



class ReadTophatFusionPre(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
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
			f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':str(self.i),'dataset':f.get_dataset_name()})# Secondary location(s)
			
			self.add_fusion(f)
	
	def parse(self):
		self.i = 0
		with open(self.filename,"r") as fh:
			for line in fh:
				self.i += 1
				self.parse_line(line)



class ReadTophatFusionPost(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		self.filename = arg_filename
		
		self.parse()
	
	def flush(self):
		self.chr_1 = False
		self.chr_2 = False
		
		self.break_1 = False
		self.break_2 = False
		
		self.seq = False
		self.insert_seq = False
	
	def parse_line_type_0(self,line):
		#self._id = line.strip()
		
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
		self.flush()
		
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
					f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':str(j),'dataset':f.get_dataset_name()})# Secondary location(s)
					
					self.add_fusion(f)
				
				i += 1



class ReadDefuse(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
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
		
		left_pos = int(line[self.parse_left_pos_column])-1
		right_pos = int(line[self.parse_right_pos_column])-1
		
		f = Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],left_pos,right_pos,line[self.parse_sequence_column],False,line[self.parse_left_strand_column],line[self.parse_right_strand_column],self.name)
		f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.get_dataset_name()})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)



class ReadChimeraScan(HighThroughputFusionDetectionExperiment):
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		self.filename = arg_filename
		self.parse()
	
	def parse_line(self,line):
		if(self.parse_header):
			self.parse_line__header(line)
		else:
			self.parse_line__fusion(line)
	
	def parse_line__header(self,line):
		line = line.strip().split("\t")
		
		self.parse_left_chr_column = line.index("#chrom5p")
		self.parse_right_chr_column = line.index("chrom3p")
		
		self.parse_left_pos_column = line.index("end5p")
		self.parse_right_pos_column = line.index("start3p")
		
		self.parse_left_strand = line.index("strand5p")
		self.parse_right_strand = line.index("strand3p")
		
		self.parse_id = line.index("chimera_cluster_id")
		
		self.parse_header = False
	
	def parse_line__fusion(self,line):
		line = line.strip().split("\t")
		
		f = Fusion(line[self.parse_left_chr_column],line[self.parse_right_chr_column],line[self.parse_left_pos_column],line[self.parse_right_pos_column],False,False,line[self.parse_left_strand],line[self.parse_right_strand],self.name)
		f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':line[self.parse_id],'dataset':f.get_dataset_name()})# Secondary location(s)
		
		self.add_fusion(f)
	
	def parse(self):
		self.parse_header = True
		
		with open(self.filename,"r") as fh:
			for line in fh:
				self.parse_line(line)





class ReadTrinityGMAP(HighThroughputFusionDetectionExperiment):
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
	
	def __init__(self,arg_filename,name):
		HighThroughputFusionDetectionExperiment.__init__(self,name,"RNA")
		
		
		self.filename = arg_filename
		self.parse()
		
	def parse(self):
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
						f.add_location({'left':[f.get_left_chromosome(True),f.get_left_break_position()],'right':[f.get_right_chromosome(True),f.get_right_break_position()],'id':uid,'dataset':f.get_dataset_name()})
						
						distance = f.get_distance()
						if(distance > 100000 or distance == -1):
							self.add_fusion(f)
					else:
						contig_chunk.append(line)
	
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

