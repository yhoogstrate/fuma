#!/usr/bin/env python



class HighThroughputFusionDetectionExperiment:
	def __init__(self,name,arg_type):
		self.name = name
		
		self.fusions_left_keys = {}
		self.fusions_right_keys = {}
		
		self.fusions_indexed_left = False
		self.fusions_indexed_right = False
		
		self.genes_overlayed_left = False
		self.genes_overlayed_right = False
		
		self.n = 0
		
		self.set_type(arg_type)
	
	def set_type(self,arg_type):
		if(arg_type in ["RNA","DNA"]):
			self.type = arg_type
		else:
			raise TypeError,"Incorrect type: "+str(arg_type)
	
	def get_type(self):
		return self.type
	
	def add_fusion(self,fusion):
		# Add left location
		left_chr = fusion.get_left_chromosome(False)
		
		if(not self.fusions_left_keys.has_key(left_chr)):
			self.fusions_left_keys[left_chr] = {}
		
		left_pos = fusion.get_left_break_position()
		
		if(not self.fusions_left_keys[left_chr].has_key(left_pos)):
			self.fusions_left_keys[left_chr][left_pos] = []
		self.fusions_left_keys[left_chr][left_pos].append(fusion)
		
		
		
		# Add right location
		right_chr = fusion.get_right_chromosome(False)
		
		if(not self.fusions_right_keys.has_key(right_chr)):
			self.fusions_right_keys[right_chr] = {}
		
		right_pos = fusion.get_right_break_position()
		
		if(not self.fusions_right_keys[right_chr].has_key(right_pos)):
			self.fusions_right_keys[right_chr][right_pos] = []
		self.fusions_right_keys[right_chr][right_pos].append(fusion)
		
		self.n += 1
	
	def count_fusions(self):
		return self.n
	
	def index_fusions_left(self):
		if(not self.converted_to_genes_left()):
			print "   - Indexing (left): "+self.name
			self.fusions_indexed_left = []
			
			for left_chr in sorted(self.fusions_left_keys.keys()):
				fusions = []
				#for location in sorted(self.fusions_left_keys[left_chr]):
				#	fusions.append(self.fusions_left_keys[left_chr][location])
				for location in sorted(self.fusions_left_keys[left_chr]):
					for fusion in self.fusions_left_keys[left_chr][location]:
						fusions.append(fusion)
				
				self.fusions_indexed_left.append({"name":left_chr,"fusions":fusions})
	
	def index_fusions_right(self):
		if(not self.converted_to_genes_right()):
			print "   - Indexing (right): "+self.name
			self.fusions_indexed_right = []
			
			for right_chr in sorted(self.fusions_right_keys.keys()):
				fusions = []
				for location in sorted(self.fusions_right_keys[right_chr]):
					fusions.append(self.fusions_right_keys[right_chr][location])
				
				self.fusions_indexed_right.append({"name":right_chr,"fusions":fusions})
	
	def find_fusions(self,fusions):
		print "Number of fusions:"
		print len(fusions)
	
	def get_fusions(self):# Make iterator object?
		"""
		Returns a SORTED list of all fusions
		"""
		if(self.fusions_indexed_left == False):
			self.index_fusions_left()
		
		return self.fusions_indexed_left
	
	def get_fusions_indexed_left(self):
		if(self.fusions_indexed_left == False):
			self.index_fusions_left()
		
		return self.fusions_indexed_left
	
	def get_fusions_indexed_right(self):
		if(self.fusions_indexed_right == False):
			self.index_fusions_right()
		
		return self.fusions_indexed_right
	
	def converted_to_genes_left(self):
		return self.genes_overlayed_left
	
	def converted_to_genes_right(self):
		return self.genes_overlayed_right
	
	def show_me(self):
		fusions = self.get_fusions_indexed_left()
		for _chr in fusions:
			for fusion in _chr["fusions"]:
				fusion.show_me()
		print "---------------------"
	
	def export_to_CG_Junctions_file(self,filename):
		fh = open(filename,"w")
		
		fh.write("#ASSEMBLY_ID	GS000007673-ASM\n")
		fh.write("#SOFTWARE_VERSION	2.0.2.20\n")
		fh.write("#GENERATED_BY	cgatools\n")
		fh.write("#GENERATED_AT	2012-Feb-16	23:05:21.058280\n")
		fh.write("#FORMAT_VERSION	2\n")
		fh.write("#GENOME_REFERENCE	NCBI	build	36\n")
		fh.write("#SAMPLE	GS00669-DNA_B02\n")
		fh.write("#TYPE	JUNCTIONS\n")
		fh.write("#DBSNP_BUILD	dbSNP	build	130\n")
		fh.write("#GENE_ANNOTATIONS	NCBI	build	36.3\n")
		fh.write("\n")
		fh.write(">Id	LeftChr	LeftPosition	LeftStrand	LeftLength	RightChr	RightPosition	RightStrand	RightLength	StrandConsistent	Interchromosomal	Distance	DiscordantMatePairAlignments	JunctionSequenceResolved	TransitionSequence	TransitionLength	LeftRepeatClassification	RightRepeatClassification	LeftGenes	RightGenes	XRef	DeletedTransposableElement	KnownUnderrepresentedRepeat	FrequencyInBaselineGenomeSet	AssembledSequence	EventId	Type	RelatedJunctions\n")
		
		fid = 1
		
		for chromosome in self.get_fusions():
			for fusion in chromosome["fusions"]:
				fh.write(str(fid)+"	")
				
				fh.write(fusion.get_left_chromosome()+"	")
				fh.write(str(fusion.get_left_break_position())+"	")
				fh.write(fusion.get_left_strand()+"	101	")
				
				fh.write(fusion.get_right_chromosome()+"	")
				fh.write(str(fusion.get_right_break_position())+"	")
				fh.write(fusion.get_right_strand()+"	101	")
				
				strand_consistent = (fusion.get_left_strand() == fusion.get_right_strand())
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
				
				if(fusion.get_transition_sequence()):
					fh.write(fusion.get_transition_sequence()+"	"+str(len(fusion.get_transition_sequence())))
				else:
					fh.write("\t")
				
				fh.write("			")
				
				fh.write(":".join(fusion.get_annotated_genes_left())+"	")
				fh.write(":".join(fusion.get_annotated_genes_right())+"	")
				
				fh.write("			1.0		"+str(fid)+"	complex	\n")
				
				fid += 1
		
		fh.close()
	
	def overlay_locations_to_genes(self,arg_dataset):
		self.overlay_left_locations_to_genes(arg_dataset)
		self.overlay_right_locations_to_genes(arg_dataset)
	
	def overlay_left_locations_to_genes(self,arg_annotations):
		if(not self.converted_to_genes_left()):
			print " - Converting LEFT breakpoints to GENES: "+self.name
			for chromosome in self.get_fusions_indexed_left():
				i = 0
				if(arg_annotations.annotations_left_indexed.has_key(chromosome["name"])):
					genes = arg_annotations.annotations_left_indexed[chromosome["name"]]
					for fusion in chromosome["fusions"]:
						current = fusion.get_annotated_genes_left()
						
						if(not current):
							current = []
							fusion.annotate_genes_left(current)
						
						for k in range(0,len(genes)-1):
							gene = genes[k]
							
							if arg_annotations.find_overlap2(fusion.get_left_break_position(),gene):
								current = fusion.get_annotated_genes_left()
								current.append(gene["name"])
								fusion.annotate_genes_left(current)
							
							k += 1
						
				
				self.genes_overlayed_left = True
		else:
			print " - Skipping"
		print "---"
	
	def overlay_right_locations_to_genes(self,arg_annotations):
		if(not self.converted_to_genes_right()):
			print " - Converting RIGHT breakpoints to GENES: "+self.name
			for chromosome in self.get_fusions_indexed_right():
				i = 0
				
				if(arg_annotations.annotations_left_indexed.has_key(chromosome["name"])):
					genes = arg_annotations.annotations_left_indexed[chromosome["name"]]
					
					for position in chromosome["fusions"]:
						for fusion in position:
							current = fusion.get_annotated_genes_right()
							
							if(not current):
								current = []
								fusion.annotate_genes_right(current)
							
							while(i < len(genes) and fusion.get_right_break_position() < genes[i]["start"]):
								i += 1
							
							k = i
							while(k < len(genes) and fusion.get_right_break_position() >= genes[k]["start"]):
								gene = genes[k]
								if arg_annotations.find_overlap2(fusion.get_right_break_position(),gene):
									current = fusion.get_annotated_genes_right()
									if(not current):
										current = []
									current.append(gene["name"])
									fusion.annotate_genes_right(current)
								
								k += 1
				
				self.genes_overlayed_right = True
		else:
			print " - *Skipping"
		print "---"
	