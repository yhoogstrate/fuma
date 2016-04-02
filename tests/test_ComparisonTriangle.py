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

import unittest,logging,sys,hashlib,os
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadDefuse
from fuma.Readers import ReadFusionMap
from fuma.Fusion import Fusion
from fuma.MergedFusion import MergedFusion
from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation
from fuma.ParseBED import ParseBED
from fuma.ComparisonTriangle import ComparisonTriangle
from fuma.CLI import CLI



def match_files_unsorted(filename_1, filename_2):
	fh1 = open(filename_1,"r")
	content1 = fh1.read().split("\n")
	fh1.close()
	
	fh2 = open(filename_2,"r")
	content2 = fh2.read().split("\n")
	fh2.close()
	
	content1.sort()
	content2.sort()
	
	if len(content1) == len(content2):
		for i in range(len(content1)):
			if content1[i] != content2[i]:
				return False
		
		return True
	else:
		return False


class TestComparisonTriangle(unittest.TestCase):
	def test_01(self):
		output_file = 'test_ComparisonTriangle.test_01.output.txt'
		validation_file = 'tests/data/test_Functional.test_01.output.txt'
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_02(self):
		output_file = 'test_ComparisonTriangle.test_02.output.txt'
		validation_file = 'tests/data/test_Functional.test_01.output.txt'
		
		args = CLI(['-m','egm','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_03(self):
		output_file = 'test_ComparisonTriangle.test_03.output.txt'
		validation_file = 'tests/data/test_Functional.test_01.output.txt'
		
		args = CLI(['-m','overlap','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	#def test_04(self):
		#"""
		#Functional test with test Edgren data (comparison to all genes on hg19)
		#"""
		
		#output_file = 'test_ComparisonTriangle.test_04.output.txt'
		#validation_file = 'tests/data/test_Functional.test_Edgren_hg19.output.list.txt'
		
		#args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o',output_file])
		
		#experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_Edgren_hg19.ChimeraScan.txt","chimerascan")
		#experiment_b = ReadDefuse("tests/data/test_Functional.test_Edgren_hg19.Defuse.txt","defuse")
		#experiment_c = ReadFusionMap("tests/data/test_Functional.test_Edgren_hg19.FusionMap.txt","fusion-map")
		#experiment_d = ReadFusionMap("tests/data/test_Functional.test_Edgren_hg19.TruePositives.txt","edgren_tp")
		
		#genes = ParseBED("tests/data/refseq_genes_hg19.bed","hg19",200000)
		
		#experiment_a.annotate_genes(genes)
		#experiment_b.annotate_genes(genes)
		#experiment_c.annotate_genes(genes)
		#experiment_d.annotate_genes(genes)
		
		#experiment_a.remove_duplicates(args)
		#experiment_b.remove_duplicates(args)
		#experiment_c.remove_duplicates(args)
		#experiment_d.remove_duplicates(args)
		
		#overlap = ComparisonTriangle(args)
		#overlap.add_experiment(experiment_a)
		#overlap.add_experiment(experiment_b)
		#overlap.add_experiment(experiment_c)
		#overlap.add_experiment(experiment_d)
		
		#overlap.overlay_fusions()
		
		#files_identical = match_files_unsorted(output_file,validation_file)
		#self.assertTrue(files_identical)
		
		#if files_identical:
			#os.remove(output_file)
	
	def test_cfbsg_01(self):
		"""
		Rewritten from test_CompareFusionsBySpanningGenes::test_01()
		"""
		output_file = 'test_ComparisonTriangle.test_cfbsg_01.output.txt'
		validation_file = None
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentB")
		
		self.assertEqual(len(experiment_a), 690)
		self.assertEqual(len(experiment_b), 690)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg18",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		
		overlap.overlay_fusions()
		
		num_lines = sum(1 for line in open(output_file,'r'))
		self.assertEqual(num_lines , 538+1)
		
		if num_lines == 538+1:
			os.remove(output_file)
	
	def test_cfbsg_02(self):
		"""
		Rewritten from test_CompareFusionsBySpanningGenes::test_02()
		"""
		
		output_file_a = 'test_ComparisonTriangle.test_cfbsg_02_a.output.txt'
		output_file_b = 'test_ComparisonTriangle.test_cfbsg_02_b.output.txt'
		
		args_a = CLI(['-m','overlap','--no-strand-specific-matching','-s','','-o',output_file_a])
		args_b = CLI(['-m','overlap',   '--strand-specific-matching','-s','','-o',output_file_b])
		
		## First test the matches if strand-specific-matching is disabled (all 4 fusions should be identical)
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02_a.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02_b.bedpe","TestExperimentB")
		
		self.assertEqual(len(experiment_a), 4)
		self.assertEqual(len(experiment_b), 4)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02.bed","hg18",200000)
		
		self.assertEqual(len(genes), 8)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		
		## @todo -> remove duplicates should be done separately
		experiment_a.remove_duplicates(args_a)
		experiment_b.remove_duplicates(args_a)
		
		# ------------------------------------ #
		# No EGM, no strand-specific-matching
		overlap = ComparisonTriangle(args_a)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		
		overlap.overlay_fusions()
		
		num_lines = sum(1 for line in open('test_ComparisonTriangle.test_cfbsg_02_a.output.txt','r'))
		self.assertEqual(num_lines , 4+1)
		
		if num_lines == 4+1:
			os.remove('test_ComparisonTriangle.test_cfbsg_02_a.output.txt')
		# ------------------------------------ #
		
		# ------------------------------------ #
		# No EGM, no strand-specific-matching
		overlap = ComparisonTriangle(args_b)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		
		overlap.overlay_fusions()
		
		num_lines = sum(1 for line in open('test_ComparisonTriangle.test_cfbsg_02_b.output.txt','r'))
		self.assertEqual(num_lines , 3+3+1+1)# 3 unmatched in A, 3 unmatched in B, 1 matched and a header line
		
		if num_lines == 1+1:
			os.remove('test_ComparisonTriangle.test_cfbsg_02_b.output.txt')
		# ------------------------------------ #
	
	def test_oc_01(self):
		output_file = 'test_ComparisonTriangle.test_oc_01.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_01.txt'
		
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_1 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment1")
		experiment_2 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment2")
		experiment_3 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment3")
		
		self.assertTrue(len(experiment_1) == 690)
		self.assertTrue(len(experiment_2) == 690)
		self.assertTrue(len(experiment_3) == 690)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_01.bed","hg18", 200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args)
		experiment_2.remove_duplicates(args)
		experiment_3.remove_duplicates(args)
		
		self.assertTrue(len(experiment_1) <= 690)
		self.assertTrue(len(experiment_2) <= 690)
		self.assertTrue(len(experiment_3) <= 690)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_02(self):
		"""
		Exp1:
		f1a: [X] -> [A,B]
		f1b: [X] -> [B,C]
		
		Exp2:
		f2a: [X] -> [A,B,C]
		
		So f2a should merge with both f1a and f1b!
		"""
		
		output_file = 'test_ComparisonTriangle.test_oc_02.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_02.txt'
		
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)#A,B
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_1","uid",True)#B,C
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_2","uid",True)#A,B,C
		fusion_4 = Fusion("chrX","chrY",10000,10000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.add_fusion(fusion_1)
		experiment_1.add_fusion(fusion_2)
		experiment_2.add_fusion(fusion_3)
		experiment_3.add_fusion(fusion_4)
		
		self.assertEqual(len(experiment_1), 2)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_02.bed","hg18", 200000)
		
		self.assertEqual(len(genes), 6)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args)
		experiment_2.remove_duplicates(args)
		experiment_3.remove_duplicates(args)
		
		self.assertEqual(len(experiment_1), 2)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_03(self):
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o','-'])
		
		fusion_1 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_1","f_ABC",True)#A,B,C
		fusion_2 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_2","f_AB",True)#A,B
		fusion_3 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_3","f_BC",True)#B,C
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_02.bed","hg18", 200000)
		
		self.assertEqual(len(genes), 6)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args)
		experiment_2.remove_duplicates(args)
		experiment_3.remove_duplicates(args)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlap = ComparisonTriangle(args)
		
		merged_fusion_12 = MergedFusion()
		merged_fusion_12.add_fusion(fusion_1)
		merged_fusion_12.add_fusion(fusion_2)
		
		merged_fusion_21 = MergedFusion()
		merged_fusion_21.add_fusion(fusion_2)
		merged_fusion_21.add_fusion(fusion_1)
		
		merged_fusion_13 = MergedFusion()
		merged_fusion_13.add_fusion(fusion_1)
		merged_fusion_13.add_fusion(fusion_3)
		
		merged_fusion_31 = MergedFusion()
		merged_fusion_31.add_fusion(fusion_3)
		merged_fusion_31.add_fusion(fusion_1)
		
		merged_fusions = [merged_fusion_12, merged_fusion_21, merged_fusion_13, merged_fusion_31]
		unique_merged_fusions = overlap.prune_duplicates(merged_fusions)
		
		self.assertEqual(len(unique_merged_fusions), 2)
	
	def test_oc_04(self):
		output_file = 'test_ComparisonTriangle.test_oc_04.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_04.txt'
		
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		gene_1 = Gene("gene_1", False)
		gene_2 = Gene("gene_2", False)
		gene_2_copy = Gene("gene_2", False)
		
		fusion_1.annotate_genes_left([gene_1])
		fusion_2.annotate_genes_left([gene_1])
		fusion_3.annotate_genes_left([gene_1])
		
		fusion_1.annotate_genes_right([gene_2])
		fusion_2.annotate_genes_right([gene_2])
		fusion_3.annotate_genes_right([gene_2_copy])
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_05(self):
		output_file = 'test_ComparisonTriangle.test_oc_05.output.txt'
		validation_file_neq = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_05.txt'
		validation_file__eq = 'tests/data/test_ComparisonTriangle.test_oc_05.output.txt'
		
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr3",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		gene_1 = Gene("gene_1", False)
		gene_2 = Gene("gene_2", False)
		
		fusion_1.annotate_genes_left([gene_1])
		fusion_2.annotate_genes_left([gene_1])
		fusion_3.annotate_genes_left([gene_1])
		
		fusion_1.annotate_genes_right([gene_2])
		fusion_2.annotate_genes_right([gene_2])
		fusion_3.annotate_genes_right([gene_2])
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		
		overlap.overlay_fusions()
		
		files_identical_1 = match_files_unsorted(output_file,validation_file__eq)
		files_identical_2 = not match_files_unsorted(output_file,validation_file_neq)
		
		"""
		This behaviour has changed since v3.0
		
		We have also found and resolved a small bug. In older versions of FuMa,
		indexing was chromosome-name based. Therefore matching two fusion genes
		only happened when they were annotated upon the same chr name. If you
		would have a fusion gene A-B (both on chr1) and fusion A-B (both on
		chr2), the old versions would consider these distinct whereas the new
		version of FuMa considers these identical.
		"""
		self.assertTrue(files_identical_1 and files_identical_2)
		
		if files_identical_1 and files_identical_2:
			os.remove(output_file)
	
	def test_oc_06(self):
		output_file = 'test_ComparisonTriangle.test_oc_06.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_06.txt'
		
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		gene_1 = Gene("gene_1", False)
		gene_2 = Gene("gene_2", False)
		gene_3 = Gene("gene_3", False)
		gene_4 = Gene("gene_4", False)
		gene_5 = Gene("gene_5", False)
		gene_6 = Gene("gene_6", False)
		
		fusion_1.annotate_genes_left([       gene_2,gene_3])
		fusion_2.annotate_genes_left([gene_1,gene_2,gene_3])
		fusion_3.annotate_genes_left([gene_4,gene_5,gene_6])
		
		fusion_1.annotate_genes_right([gene_4,gene_5,gene_6])
		fusion_2.annotate_genes_right([gene_4,gene_5       ])
		fusion_3.annotate_genes_right([gene_1,gene_2,gene_3])
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_07(self):
		"""
		Experiment1:
		f1: [X] -> [A,B]
		
		Experiment2:
		f2: [X] -> [B,C]
		
		Experiment3:
		f1: [X] -> [A,B,C]
		
		
		A,B is a subset of A,B,C
		B,C is a subset of A,B,C
		A,B is NOT a subset of B,C
		
		This means that the final overlap between all 3 experiments must be 0 independent of the order of matching
		->
		n overlap = 0
		"""
		output_file_1 = 'test_ComparisonTriangle.test_oc_07_1.output.txt'
		output_file_2 = 'test_ComparisonTriangle.test_oc_07_2.output.txt'
		output_file_3 = 'test_ComparisonTriangle.test_oc_07_3.output.txt'
		output_file_4 = 'test_ComparisonTriangle.test_oc_07_4.output.txt'
		output_file_5 = 'test_ComparisonTriangle.test_oc_07_5.output.txt'
		output_file_6 = 'test_ComparisonTriangle.test_oc_07_6.output.txt'
		validation_file_1 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_1.txt'
		validation_file_2 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_2.txt'
		validation_file_3 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_3.txt'
		validation_file_4 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_4.txt'
		validation_file_5 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_5.txt'
		validation_file_6 = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_07_6.txt'
		
		args_1 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_1])
		args_2 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_2])
		args_3 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_3])
		args_4 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_4])
		args_5 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_5])
		args_6 = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s','','-o',output_file_6])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_07.bed","hg18", 200000)
		
		self.assertEqual(len(genes), 4)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args_1)
		experiment_2.remove_duplicates(args_1)
		experiment_3.remove_duplicates(args_1)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		
		overlap_1 = ComparisonTriangle(args_1)
		overlap_1.add_experiment(experiment_1)
		overlap_1.add_experiment(experiment_2)
		overlap_1.add_experiment(experiment_3)
		
		overlap_1.overlay_fusions()
		files_identical = match_files_unsorted(output_file_1,validation_file_1)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_1)
		
		
		overlap_2 = ComparisonTriangle(args_2)
		overlap_2.add_experiment(experiment_1)
		overlap_2.add_experiment(experiment_3)
		overlap_2.add_experiment(experiment_2)
		
		overlap_2.overlay_fusions()
		files_identical = match_files_unsorted(output_file_2,validation_file_2)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_2)
		
		
		overlap_3 = ComparisonTriangle(args_3)
		overlap_3.add_experiment(experiment_2)
		overlap_3.add_experiment(experiment_1)
		overlap_3.add_experiment(experiment_3)
		
		overlap_3.overlay_fusions()
		files_identical = match_files_unsorted(output_file_3,validation_file_3)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_3)
		
		
		overlap_4 = ComparisonTriangle(args_4)
		overlap_4.add_experiment(experiment_2)
		overlap_4.add_experiment(experiment_3)
		overlap_4.add_experiment(experiment_1)
		
		overlap_4.overlay_fusions()
		files_identical = match_files_unsorted(output_file_4,validation_file_4)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_4)
		
		
		overlap_5 = ComparisonTriangle(args_5)
		overlap_5.add_experiment(experiment_3)
		overlap_5.add_experiment(experiment_1)
		overlap_5.add_experiment(experiment_2)
		
		overlap_5.overlay_fusions()
		files_identical = match_files_unsorted(output_file_5,validation_file_5)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_5)

		
		overlap_6 = ComparisonTriangle(args_6)
		overlap_6.add_experiment(experiment_3)
		overlap_6.add_experiment(experiment_2)
		overlap_6.add_experiment(experiment_1)
		
		overlap_6.overlay_fusions()
		files_identical = match_files_unsorted(output_file_6,validation_file_6)
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_file_6)
	
	def test_oc_08(self):
		"""
		This tests whether intergenic fusions are taken into account.
		Both input files contain one intergenic fusion, from gene1 to gene1.
		
		The genomic location of the gene is chr1:15000-16000 and the
		breakpoints of the fusion are: chr1:15250 and chr1:15750.
		
		In the first tests, it checks whether this fusion gene (chr1:15250-chr1:15750)
		with on both genes "gene1" annotated, picked up by two experiments,
		is indeed considered identical.
		
		In the second test, also the strands are taken into account. Because
		in the first experiment, the strands are both positive and in the
		second experiment the strands are both negative, FuMa should not
		match them and the size of the corresponding matched dataset should
		be 0.
		"""
		output_file_a = 'test_ComparisonTriangle.test_oc_08_a.output.txt'
		output_file_b = 'test_ComparisonTriangle.test_oc_08_b.output.txt'
		validation_file_a = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_08_a.txt'
		validation_file_b = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_08_b.txt'
		
		args_a = CLI(['-m','subset','--no-strand-specific-matching','-f','list','-s','','-o',output_file_a])
		args_b = CLI(['-m','subset',   '--strand-specific-matching','-f','list','-s','','-o',output_file_b])
		
		experiment_1 = ReadFusionMap("tests/data/test_OverlapComplex.TestOverlapComplex.test_08_ff.FusionMap.txt","TestExperiment1")
		experiment_2 = ReadFusionMap("tests/data/test_OverlapComplex.TestOverlapComplex.test_08_rr.FusionMap.txt","TestExperiment2")
		
		self.assertTrue(len(experiment_1) == 1)
		self.assertTrue(len(experiment_2) == 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_08.bed","hg19", 200000)
		
		self.assertEqual(len(genes), 1)# 1 gene; intergenic fusion
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args_a)
		experiment_2.remove_duplicates(args_a)
		
		self.assertTrue(len(experiment_1) == 1)
		self.assertTrue(len(experiment_2) == 1)
		
		# Matching, do not take strand-specific-matching into account
		overlap = ComparisonTriangle(args_a)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file_a,validation_file_a)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file_a)
		
		# Matching, do take strand-specific-matching into account.
		overlap = ComparisonTriangle(args_b)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file_b,validation_file_b)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file_b)
	
	def test_oc_09(self):
		"""
		Tests whether the overlap() matching function is implemented correctly 

Following exammple

    1200 1400 1600 1800
      :    :   :    :
      f1   f2  f3   f4
      |    |   |    |
[---A1--]  |   |    |
     [---A2--] |    |
         [---A3--]  |
             [---A4--]
    [A5]     [-----A5---]
                   [---A6--]


f1=[--A1--],[--A2--],                  [--A5--]
f2=         [--A2--],[--A3--]
f3=                  [--A3--],[--A4--],[--A5--]
f4=                           [--A4--],[--A5--],[--A6--]

exp1=[f1]
exp1=[f2]
exp1=[f3]
exp1=[f4]

exp1,exp2 = ([--A2--])
exp1,exp3 = ([--A5--])
exp1,exp4 = ([--A5--])
exp2,exp3 = ([--A3--])
exp2,exp4 = none
exp3,exp4 = ([--A4--],[--A5--])


(exp1,exp2),exp3 = ([--A2--])          , ([--A3--],[--A4--],[--A5--]) = none
(exp1,exp3),exp2 = ([--A5--])          , ([--A2--],[--A3--])          = none
(exp2,exp3),exp1 = ([--A3--])          , ([--A1--],[--A2--],[--A5--]) = none

(exp1,exp2),exp4 = ([--A2--])          , ([--A4--],[--A5--],[--A6--]) = none
(exp1,exp4),exp2 = ([--A5--])          , ([--A2--],[--A3--])          = none
(exp2,exp4),exp1 = none                , ([--A1--],[--A2--],[--A5--]) = none

(exp1,exp3),exp4 = ([--A5--])          , ([--A4--],[--A5--],[--A6--]) = ([--A5--])
(exp1,exp4),exp3 = ([--A5--])          , ([--A1--],[--A2--],[--A5--]) = ([--A5--])
(exp3,exp4),exp1 = ([--A4--],[--A5--]) , ([--A1--],[--A2--],[--A5--]) = ([--A5--])

(exp2,exp3),exp4 = ([--A3--])          , ([--A4--],[--A5--],[--A6--]) = none
(exp2,exp4),exp3 = none                , ([--A3--],[--A4--],[--A5--]) = none
(exp3,exp4),exp2 = ([--A4--],[--A5--]) , ([--A2--],[--A3--])          = none

unique fusions
(exp1,exp2)
(exp2,exp3)
(exp1,exp3,exp4)
		"""
		output_files = []
		validation_files = []
		args = []
		for i in range(1,19+1):
			str_i = str(i)
			while len(str_i) < 2:
				str_i = "0"+str_i
			
			output = "test_ComparisonTriangle.test_oc_09_"+str_i+".output.txt"
			
			output_files.append(output)
			validation_files.append("tests/data/test_OverlapComplex.TestOverlapComplex.test_09_"+str_i+".output.txt")
			args.append(CLI(['-m','overlap','-f','list','--strand-specific-matching','-s','','-o',output]))
		
		genes = GeneAnnotation("hg19")
		gene_A1 = Gene("[--A1--]", False)
		gene_A2 = Gene("[--A2--]", False)
		gene_A3 = Gene("[--A3--]", False)
		gene_A4 = Gene("[--A4--]", False)
		gene_A5 = Gene("[--A5--]", False)
		gene_A6 = Gene("[--A6--]", False)
		gene_XX = Gene("X", False)
		
		genes.add_annotation(gene_A1,"1",10000,13000)
		genes.add_annotation(gene_A2,"1",11500,14500)
		genes.add_annotation(gene_A3,"1",13000,17000)
		genes.add_annotation(gene_A4,"1",15000,18500)
		genes.add_annotation(gene_A5,"1",15000,19000)
		genes.add_annotation(gene_A5,"1",11500,12500)# Add twice
		genes.add_annotation(gene_A6,"1",17000,19000)
		genes.add_annotation(gene_XX,"X",14000,16000)
		
		fusion_1 = Fusion("chr1","chrX",12000,15000,"+","+","Experiment_1","1",True)
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chrX",14000,15000,"+","+","Experiment_2","2",True)
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		fusion_3 = Fusion("chr1","chrX",16000,15000,"+","+","Experiment_3","3",True)
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		experiment_3.add_fusion(fusion_3)
		experiment_3.annotate_genes(genes)
		
		fusion_4 = Fusion("chr1","chrX",18000,15000,"+","+","Experiment_4","4",True)
		experiment_4 = FusionDetectionExperiment("Experiment_4")
		experiment_4.add_fusion(fusion_4)
		experiment_4.annotate_genes(genes)
		
		i = 0
		
		#(b1,b2) = (A2)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3) = (A5)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4) = (A2)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3) = (A3)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4) = (A4,A5)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2),b3 = (A2),(A3,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3),b2 = (A5),(A2,A3)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3),b1 = (A3),(A1,A2,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2),b4 = (A2),(A4,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4),b2 = (A5),(A2,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4),b1 = (A4),(A1,A2,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3),b4 = (A3),(A4,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4),b3 = (A5),(A3,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4),b1 = (A4),(A1,A3,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3),b4 = (A3),(A4,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4),b3 = (A5),(A3,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4),b2 = (A4),(A2,A3,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2,b3,b4) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
	
	def test_oc_10(self):
		"""
		Tests whether the overlap() matching function is implemented correctly 

Following exammple

    1200 1400 1600 1800
      :    :   :    :
      f1   f2  f3   f4
      |    |   |    |
[---A1--]  |   |    |
     [---A2--] |    |
         [---A3--]  |
          [---A4--] |
            [---A5---]
                   [---A6--]


f1=[--A1--],[--A2--]
f2=         [--A2--],[--A3--],[--A4--]
f3=                  [--A3--],[--A4--],[--A5--]
f4=                                    [--A5--],[--A6--]

exp1 = [f1] = ([--A1--],[--A2--])
exp1 = [f2] = ([--A2--],[--A3--],[--A4--])
exp1 = [f3] = ([--A3--],[--A4--],[--A5--])
exp1 = [f4] = ([--A5--],[--A6--])

exp1,exp2 = [--A2--]
exp1,exp3 = none
exp1,exp4 = none
exp2,exp3 = [--A3--],[--A4--]
exp2,exp4 = none
exp3,exp4 = [--A5--]


# 1,2,3
(exp1,exp2),exp3 = ([--A2--])          , ([--A3--],[--A4--],[--A5--]) = none
(exp1,exp3),exp2 = (none)              , ([--A2--],[--A3--],[--A4--]) = none
(exp2,exp3),exp1 = ([--A3--],[--A4--]) , ([--A1--],[--A2--])          = none

# 1,2,4
(exp1,exp2),exp4 = ([--A2--])          , ([--A5--],[--A6--])          = none
(exp1,exp4),exp2 = (none)              , ([--A2--],[--A3--],[--A4--]) = none
(exp2,exp4),exp1 = (none)              , ([--A1--],[--A2--])          = none

# 1,3,4
(exp1,exp3),exp4 = (none)              , ([--A5--],[--A6--])          = none
(exp1,exp4),exp3 = (none)              , ([--A3--],[--A4--],[--A5--]) = none
(exp3,exp4),exp1 = ([--A5--])          , ([--A1--],[--A2--])          = none

# 2,3,4
(exp2,exp3),exp4 = ([--A3--],[--A4--]) , ([--A5--],[--A6--])          = none
(exp2,exp4),exp3 = (none)              , ([--A3--],[--A4--],[--A5--]) = none
(exp3,exp4),exp2 = ([--A5--])          , ([--A2--],[--A3--],[--A4--]) = none

unique fusions
(exp1,exp2): [--A2--]
(exp2,exp3): [--A3--],[--A4--]
(exp3,exp4): [--A5--]
		"""
		
		output_files = []
		validation_files = []
		args = []
		for i in range(1,19+1):
			str_i = str(i)
			while len(str_i) < 2:
				str_i = "0"+str_i
			
			output = "test_ComparisonTriangle.test_oc_10_"+str_i+".output.txt"
			
			output_files.append(output)
			validation_files.append("tests/data/test_OverlapComplex.TestOverlapComplex.test_10_"+str_i+".output.txt")
			args.append(CLI(['-m','overlap','-f','list','--strand-specific-matching','-s','','-o',output]))
		
		genes = GeneAnnotation("hg19")
		gene_A1 = Gene("[--A1--]", False)
		gene_A2 = Gene("[--A2--]", False)
		gene_A3 = Gene("[--A3--]", False)
		gene_A4 = Gene("[--A4--]", False)
		gene_A5 = Gene("[--A5--]", False)
		gene_A6 = Gene("[--A6--]", False)
		gene_XX = Gene("X", False)
		
		genes.add_annotation(gene_A1,"1",10000,13000)
		genes.add_annotation(gene_A2,"1",11500,14500)
		genes.add_annotation(gene_A3,"1",13000,17000)
		genes.add_annotation(gene_A4,"1",13500,17500)
		genes.add_annotation(gene_A5,"1",15000,18500)
		genes.add_annotation(gene_A6,"1",17500,19000)
		genes.add_annotation(gene_XX,"X",14000,16000)
		
		fusion_1 = Fusion("chr1","chrX",12000,15000,"+","+","Experiment_1","1",True)
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chrX",14000,15000,"+","+","Experiment_2","2",True)
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		fusion_3 = Fusion("chr1","chrX",16000,15000,"+","+","Experiment_3","3",True)
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		experiment_3.add_fusion(fusion_3)
		experiment_3.annotate_genes(genes)
		
		fusion_4 = Fusion("chr1","chrX",18000,15000,"+","+","Experiment_4","4",True)
		experiment_4 = FusionDetectionExperiment("Experiment_4")
		experiment_4.add_fusion(fusion_4)
		experiment_4.annotate_genes(genes)
		
		i = 0
		
		#(b1,b2) = (A2)
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3) = [--A3--],[--A4--]
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4) = [--A5--]
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2),b3 = (A2),(A3,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3),b2 = (A5),(A2,A3)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3),b1 = (A3),(A1,A2,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2),b4 = (A2),(A4,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4),b2 = (A5),(A2,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4),b1 = (A4),(A1,A2,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b3),b4 = (A3),(A4,A4,A5)    = none
		print output_files[i]
		print validation_files[i]

		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b4),b3 = (A5),(A3,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4),b1 = (A4),(A1,A3,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_1)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b3),b4 = (A3),(A4,A4,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b2,b4),b3 = (A5),(A3,A4)       = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b3,b4),b2 = (A4),(A2,A3,A5)    = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
		i += 1
		
		#(b1,b2,b3,b4) = none
		overlap = ComparisonTriangle(args[i])
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_files[i],validation_files[i])
		self.assertTrue(files_identical)
		if files_identical:
			os.remove(output_files[i])
	
	def test_oc_11(self):
		"""
		Tests whether the overlap() matching function is implemented correctly 

Following exammple

  1200 1400 1600 1800      5000
   :    :      :    :       :
   f1   f2     f3   f4      f5
   |    |      |    |       |
[-A1-]  |      |  [-A1-]    |
[----A2----]   |            |
        |   [----A3----]    |
        |   [----A4----]    |
      [----A5----]          |
                          [ A6 ]

f1=[--A1--],[--A2--]
f2=[        [--A2--],                  [--A5--]
f3=                  [--A3--],[--A4--]
f4=[--A1--]          [--A3--],[--A4--]
f5=                                            [--A6--]
		"""
		output_file = 'test_ComparisonTriangle.test_oc_11.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_11.output.txt'
		
		args = CLI(['-m','overlap','-f','list','--strand-specific-matching','-s','','-o',output_file])
		
		genes = GeneAnnotation("hg19")
		gene_A1 = Gene("[--A1--]", False)
		gene_A2 = Gene("[--A2--]", False)
		gene_A3 = Gene("[--A3--]", False)
		gene_A4 = Gene("[--A4--]", False)
		gene_A5 = Gene("[--A5--]", False)
		gene_A6 = Gene("[--A6--]", False)
		gene_XX = Gene("X", False)
		
		genes.add_annotation(gene_A1,"1",10000,13000)
		genes.add_annotation(gene_A1,"1",17000,19000)
		genes.add_annotation(gene_A2,"1",10000,14800)
		genes.add_annotation(gene_A3,"1",15200,20000)
		genes.add_annotation(gene_A4,"1",15200,20000)
		genes.add_annotation(gene_A5,"1",13000,17000)
		genes.add_annotation(gene_A6,"1",49000,51000)
		genes.add_annotation(gene_XX,"X",14000,16000)
		
		
		fusion_1 = Fusion("chr1","chrX",12000,15000,"+","+","Experiment_1","1",True)
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chrX",14000,15000,"+","+","Experiment_2","2",True)
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		fusion_3 = Fusion("chr1","chrX",16000,15000,"+","+","Experiment_3","3",True)
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		experiment_3.add_fusion(fusion_3)
		experiment_3.annotate_genes(genes)
		
		fusion_4 = Fusion("chr1","chrX",18000,15000,"+","+","Experiment_4","4",True)
		experiment_4 = FusionDetectionExperiment("Experiment_4")
		experiment_4.add_fusion(fusion_4)
		experiment_4.annotate_genes(genes)
		
		fusion_5 = Fusion("chr1","chrX",50000,15000,"+","+","Experiment_5","5",True)
		experiment_5 = FusionDetectionExperiment("Experiment_5")
		experiment_5.add_fusion(fusion_5)
		experiment_5.annotate_genes(genes)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.add_experiment(experiment_4)
		overlap.add_experiment(experiment_5)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_12(self):
		"""
		Tests whether the overlap() matching function is implemented correctly 
- This is the two examples given in the github manual -
		"""
		output_file = 'test_ComparisonTriangle.test_oc_12.output.txt'
		validation_file = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_12.output.txt'
		
		args = CLI(['-m','overlap','-f','list','-s','','-o',output_file])
		
		gene_green  = Gene("GREEN", False)
		gene_blue   = Gene("BLUE", False)
		gene_yellow = Gene("YELLOW", False)
		gene_purple = Gene("PURPLE", False)
		gene_XX     = Gene("X", False)
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_blue,  "1",12000,14000)
		genes.add_annotation(gene_green, "1",13000,14000)
		genes.add_annotation(gene_yellow,"1",16000,18000)
		genes.add_annotation(gene_purple,"1",12000,13000)
		
		fusion_1 = Fusion("chr1","chr1",12500,17000,"+","+","Experiment_1","uid",True)
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		
		fusion_2 = Fusion("chr1","chr1",13500,17000,"+","+","Experiment_2","uid",True)
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file,validation_file)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file)
	
	def test_oc_13(self):
		""" 
	f1 = A,B,C
	f2 = A,B
	f3 = B,C
	
     f2  f1 f3
     |   |  |
[ -- A -- ]
   [ -- B -- ]
      [ -- C -- ]
	
	What will be the outcome? order dependent?
	
	-> both for 'overlap' and 'subset' based matching
		"""
		output_file_a = 'test_ComparisonTriangle.test_oc_13_subset.output.txt'
		output_file_b = 'test_ComparisonTriangle.test_oc_13_overlap.output.txt'
		validation_file_a = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_13_subset.output.txt'
		validation_file_b = 'tests/data/test_OverlapComplex.TestOverlapComplex.test_13_overlap.output.txt'
		
		args_subset  = CLI(['-m', 'subset','-f','list','-s','','-o',output_file_a])
		args_overlap = CLI(['-m','overlap','-f','list','-s','','-o',output_file_b])
		
		gene_A = Gene("A", False)
		gene_B = Gene("B", False)
		gene_C = Gene("C", False)
		gene_X = Gene("X", False)
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",12000,16000)
		genes.add_annotation(gene_B,"1",13000,17000)
		genes.add_annotation(gene_C,"1",14000,18000)
		genes.add_annotation(gene_X,"X",10000,20000)
		
		fusion_1 = Fusion("chr1","chrX",15000,15000,"+","+","Experiment_1","1",True)
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chrX",13500,15000,"+","+","Experiment_2","2",True)
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		fusion_3 = Fusion("chr1","chrX",16500,15000,"+","+","Experiment_3","3",True)
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		experiment_3.add_fusion(fusion_3)
		experiment_3.annotate_genes(genes)
		
		self.assertEqual(len(fusion_1.annotated_genes_left) , 3)
		self.assertEqual(len(fusion_2.annotated_genes_left) , 2)
		self.assertEqual(len(fusion_3.annotated_genes_left) , 2)
		
		#---------------------------------------------------------------#
		# Matching, do not take strand-specific-matching into account
		overlap = ComparisonTriangle(args_subset)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file_a,validation_file_a)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file_a)
		#---------------------------------------------------------------#
		
		#---------------------------------------------------------------#
		# Matching, do take strand-specific-matching into account.
		overlap = ComparisonTriangle(args_overlap)
		overlap.add_experiment(experiment_1)
		overlap.add_experiment(experiment_2)
		overlap.add_experiment(experiment_3)
		overlap.overlay_fusions()
		
		files_identical = match_files_unsorted(output_file_b,validation_file_b)
		self.assertTrue(files_identical)
		
		if files_identical:
			os.remove(output_file_b)
		#---------------------------------------------------------------#


def main():
	unittest.main()

if __name__ == '__main__':
	main()
