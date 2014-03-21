#!/usr/bin/env python



import sys,os,os.path,argparse



sys.path.append("fuma")
sys.path.append("fuma/fuma")



from ParseBED import ParseBED
from OverlayFusions import OverlayFusions

from Readers import *

from CompareFusionsGTFOverlay import CompareFusionsGTFOverlay



if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-a","--add-gene-annotation",help="alias:filename  * file in BED format",nargs="*")
	
	parser.add_argument("-s","--add-sample",help="alias:type:filename",nargs="+")
	parser.add_argument("-l","--link-sample-to-annotation",help="sample_alias:annotation_alias",nargs="*")
	parser.add_argument("-o","--output",help="output filename; '-' for stdout",default="overlap/")
	
	args = parser.parse_args()
	
	gene_annotations = {}
	if(args.add_gene_annotation):
		for gene_annotation in args.add_gene_annotation:
			gene_annotation = gene_annotation.split(":",1)
			gene_annotations[gene_annotation[0]] = ParseBED(gene_annotation[1],gene_annotation[0])
	
	samples = {}
	sample_names = []
	for sample in args.add_sample:
		
		sample = sample.split(":",2)
		
		if(sample[0] in sample_names):
			raise Exception("non-unique sample alias: "+sample[0])
		elif(sample[0].find("~") > -1):
			raise Exception("a sample alias may not include the '~' char: "+sample[0])
		else:
			sample_names.append(sample[0])
			sample[1] = sample[1].lower().replace("-","").replace("_","").replace(" ","")
			
			if(sample[1] in ["cg","completegenomics"]):
				samples[sample[0]] = ReadCGhighConfidenceJunctionsBeta(sample[2],sample[0])
			elif(sample[1] == "chimerascan"):
				samples[sample[0]] = ReadChimeraScan(sample[2],sample[0])
			elif(sample[1] == "defuse"):
				samples[sample[0]] = ReadDefuse(sample[2],sample[0])
			elif(sample[1] in ["illuminahiseq","illuminahiseqvcf"]):
				samples[sample[0]] = ReadIlluminaHiSeqVCF(sample[2],sample[0])
			elif(sample[1] == "tophatfusionpost"):
				samples[sample[0]] = ReadTophatFusionPost(sample[2],sample[0])
			elif(sample[1] == "tophatfusionpre"):
				samples[sample[0]] = ReadTophatFusionPre(sample[2],sample[0])
			elif(sample[1] == "trinitygmap"):
				samples[sample[0]] = ReadTrinityGMAP(sample[2],sample[0])
			else:
				raise Exception("unsupported sample type: "+sample[1])
	
	if(args.link_sample_to_annotation):
		tmp = CompareFusionsGTFOverlay(False,False)
		
		for link in args.link_sample_to_annotation:
			link = link.split(":",1)
			
			if(samples.has_key(link[0])):
				sample = samples[link[0]]
			else:
				raise Exception("unknown sample: "+link[0])
			
			if(gene_annotations.has_key(link[1])):
				gene_annotation = gene_annotations[link[1]]
			else:
				raise Exception("unknown annotation: "+link[1])
			
			sample.overlay_left_locations_to_genes(gene_annotation)
			sample.overlay_right_locations_to_genes(gene_annotation)
			
			samples[link[0]] = tmp.remove_duplicates(sample)
	
	o = OverlayFusions()
	
	for sample_name in sample_names:
		o.add_dataset(samples[sample_name])
	
	o.overlay_fusions()
	o.export3(args.output)


