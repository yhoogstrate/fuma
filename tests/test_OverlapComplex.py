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

import unittest

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadFusionMap
from fuma.ParseBED import ParseBED
from fuma.OverlapComplex import OverlapComplex
from fuma.Fusion import Fusion
from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation

class TestOverlapComplex(unittest.TestCase):
	# The test indicated a minor bug. The unittest aborts on this test
	# case. A solution is given at the following url:
	# <http://stackoverflow.com/questions/4732827/continuing-in-pythons-unittest-when-an-assertion-fails>
	
	#def setUp(self):
	#	self.verificationErrors = []
	
	#def tearDown(self):
	#	self.assertEqual([], self.verificationErrors)
	
	def test_01(self):
		experiment_1 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment1")
		experiment_2 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment2")
		experiment_3 = ReadChimeraScanAbsoluteBEDPE("tests/data/test_OverlapComplex.TestOverlapComplex.test_01.bedpe","TestExperiment3")
		
		self.assertTrue(len(experiment_1) == 690)
		self.assertTrue(len(experiment_2) == 690)
		self.assertTrue(len(experiment_3) == 690)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_01.bed","hg18")
		
		self.assertEqual(len(genes), 47790)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates("by-gene-names")
		experiment_2.remove_duplicates("by-gene-names")
		experiment_3.remove_duplicates("by-gene-names")
		
		self.assertTrue(len(experiment_1) <= 690)
		self.assertTrue(len(experiment_2) <= 690)
		self.assertTrue(len(experiment_3) <= 690)
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,False,"summary")
		
		self.assertLessEqual(len(overlap), len(experiment_1))
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 538)
	
	def test_02(self):
		"""
		Exp1:
		f1: [X] -> [A,B]
		f2: [X] -> [B,C]
		
		Exp2:
		f1: [X] -> [A,B,C]
		
		Expected Exp1+2:
		f1: [X] -> [B] << subset of (A,B)(B,C)
		
		n overlap = 1 (overlap is measured as the merged fusions - only
		one merged fusion can account for all of them)
		
		                               Experiment_1 (2)              Experiment_2 (1)         Experiment_3 (1)
		Experiment_1 (2)                                  2/2 (100.0%) : 2/1 (200.0%)  0/2 (0.0%) : 0/1 (0.0%)
		Experiment_2 (1)    2/1 (200.0%) : 2/2 (100.0%)                                0/1 (0.0%) : 0/1 (0.0%)
		Experiment_3 (1)        0/1 (0.0%) : 0/2 (0.0%)       0/1 (0.0%) : 0/1 (0.0%)                         

		
		                                Experiment_1 (2)  Experiment_2 (1)         Experiment_3 (1)
		Experiment_1 & Experiment_2 (2)                                     0/2 (0.0%) : 0/1 (0.0%)
		Experiment_1 & Experiment_3 (0)                     0 : 0/1 (0.0%)  
		Experiment_2 & Experiment_3 (0)   0 : 0/2 (0.0%)                      
		"""
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")#A,B,C
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_1")#B,C
		fusion_3 = Fusion("chrX","chr2",15000,70000,None,None,"+","+","Experiment_2")#A,B
		fusion_4 = Fusion("chrX","chrY",10000,10000,None,None,"+","+","Experiment_3")
		
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		
		experiment_1.add_fusion(fusion_1)
		experiment_1.add_fusion(fusion_2)
		experiment_2.add_fusion(fusion_3)
		experiment_3.add_fusion(fusion_4)
		
		self.assertEqual(len(experiment_1), 2)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02.bed","hg18")
		
		self.assertEqual(len(genes), 6)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates("by-gene-names")
		experiment_2.remove_duplicates("by-gene-names")
		experiment_3.remove_duplicates("by-gene-names")
		
		self.assertEqual(len(experiment_1), 2)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,False,"summary")
		
		#overlapping_complex.export_summary("-")
		
		self.assertEqual(overlapping_complex.matches_total['1.2'],2)
		self.assertEqual(overlapping_complex.matches_total['1.2.3'],0)
	
	def test_04(self):
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_2")
		fusion_3 = Fusion("chrX","chr2",15000,70000,None,None,"+","+","Experiment_3")
		
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
		gene_1 = Gene("gene_1")
		gene_2 = Gene("gene_2")
		gene_2_copy = Gene("gene_2")
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,False,"summary")
		
		#overlapping_complex.export_summary("-")
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 1)
	
	def test_05(self):
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_2")
		fusion_3 = Fusion("chrX","chr3",15000,70000,None,None,"+","+","Experiment_3")
		
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
		gene_1 = Gene("gene_1")
		gene_2 = Gene("gene_2")
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,False,"summary")
		
		#overlapping_complex.export_summary("-")
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == 1)
		self.assertTrue(overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 0)
	
	def test_06(self):
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_2")
		fusion_3 = Fusion("chrX","chr2",15000,70000,None,None,"+","+","Experiment_3")
		
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
		gene_1 = Gene("gene_1")
		gene_2 = Gene("gene_2")
		gene_3 = Gene("gene_3")
		gene_4 = Gene("gene_4")
		gene_5 = Gene("gene_5")
		gene_6 = Gene("gene_6")
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,False,"summary")
		
		#overlapping_complex.export_summary("-")
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == 1)
		self.assertTrue(overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 0)
	
	def test_07(self):
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
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_2")
		fusion_3 = Fusion("chrX","chr2",15000,70000,None,None,"+","+","Experiment_3")
		
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		
		experiment_1.add_fusion(fusion_1)
		experiment_2.add_fusion(fusion_2)
		experiment_3.add_fusion(fusion_3)
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_07.bed","hg18")
		
		self.assertEqual(len(genes), 4)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates("by-gene-names")
		experiment_2.remove_duplicates("by-gene-names")
		experiment_3.remove_duplicates("by-gene-names")
		
		self.assertEqual(len(experiment_1), 1)
		self.assertEqual(len(experiment_2), 1)
		self.assertEqual(len(experiment_3), 1)
		
		overlapping_complex_1 = OverlapComplex()
		overlapping_complex_1.add_experiment(experiment_1)
		overlapping_complex_1.add_experiment(experiment_2)
		overlapping_complex_1.add_experiment(experiment_3)
		
		overlapping_complex_2 = OverlapComplex()
		overlapping_complex_2.add_experiment(experiment_1)
		overlapping_complex_2.add_experiment(experiment_3)
		overlapping_complex_2.add_experiment(experiment_2)
		
		overlapping_complex_3 = OverlapComplex()
		overlapping_complex_3.add_experiment(experiment_2)
		overlapping_complex_3.add_experiment(experiment_1)
		overlapping_complex_3.add_experiment(experiment_3)
		
		overlapping_complex_4 = OverlapComplex()
		overlapping_complex_4.add_experiment(experiment_2)
		overlapping_complex_4.add_experiment(experiment_3)
		overlapping_complex_4.add_experiment(experiment_1)
		
		overlapping_complex_5 = OverlapComplex()
		overlapping_complex_5.add_experiment(experiment_3)
		overlapping_complex_5.add_experiment(experiment_1)
		overlapping_complex_5.add_experiment(experiment_2)
		
		overlapping_complex_6 = OverlapComplex()
		overlapping_complex_6.add_experiment(experiment_3)
		overlapping_complex_6.add_experiment(experiment_2)
		overlapping_complex_6.add_experiment(experiment_1)
		
		overlap_1 = overlapping_complex_1.overlay_fusions(True,False,"summary")
		overlap_2 = overlapping_complex_2.overlay_fusions(True,False,"summary")
		overlap_3 = overlapping_complex_3.overlay_fusions(True,False,"summary")
		overlap_4 = overlapping_complex_4.overlay_fusions(True,False,"summary")
		overlap_5 = overlapping_complex_5.overlay_fusions(True,False,"summary")
		overlap_6 = overlapping_complex_6.overlay_fusions(True,False,"summary")
		
		self.assertEqual(overlapping_complex_1.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_2.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_3.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_4.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_5.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_6.matches_total['1.2.3'],0)
	
	def test_08(self):
		"""This tests whether intergenic fusions are taken into account.
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
		
		experiment_1 = ReadFusionMap("tests/data/test_OverlapComplex.TestOverlapComplex.test_08_ff.FusionMap.txt","TestExperiment1")
		experiment_2 = ReadFusionMap("tests/data/test_OverlapComplex.TestOverlapComplex.test_08_rr.FusionMap.txt","TestExperiment2")
		
		self.assertTrue(len(experiment_1) == 1)
		self.assertTrue(len(experiment_2) == 1)
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_08.bed","hg19")
		
		self.assertEqual(len(genes), 1)# 1 gene; intergenic fusion
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		
		experiment_1.remove_duplicates("by-gene-names")
		experiment_2.remove_duplicates("by-gene-names")
		
		self.assertTrue(len(experiment_1) == 1)
		self.assertTrue(len(experiment_2) == 1)
		
		# Matching, do not take strand-specific-matching into account
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=False)
		self.assertTrue(len(overlap[0]) == 1)
		
		# Matching, do take strand-specific-matching into account.
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
	
	def test_09(self):
		""" Tests whether the overlap() matching function is implemented correctly 

Following exammple

      f1   f2  f3   f4
      |    |   |    |
[---A1--]  |   |    |
     [---A2--] |    |
         [---A3--]  |
             [---A4--]
    [A5]     [-----A5---]
                   [---A6--]

b1=A1,A2,      A5
b2=   A2,A3
b3=      A3,A4,A5
b4=         A4,A5,A6

b1,b2 = (A2)
b1,b3 = (A5)
b1,b4 = (A5)
b2,b3 = (A3)
b2,b4 = none
b3,b4 = (A4,A5)


(b1,b2),b3 = (A2),(A3,A4,A5)    = none
(b1,b3),b2 = (A5),(A2,A3)       = none
(b2,b3),b1 = (A3),(A1,A2,A5)    = none

(b1,b2),b4 = (A2),(A4,A5,A6)    = none
(b1,b4),b2 = (A5),(A2,A3)       = none
(b2,b4),b1 = none,(A1,A2,A5)    = none

(b1,b3),b4 = (A5),(A4,A5,A6)    = (A5)
(b1,b4),b3 = (A5),(A1,A2,A5)    = (A5)
(b3,b4),b1 = (A4,A5),(A1,A2,A5) = (A5)

(b2,b3),b4 = (A3),(A4,A5,A6)    = none
(b2,b4),b3 = none,(A3,A4,A5)    = none
(b3,b4),b2 = (A4,A5),(A2,A3)    = none

unique fusions
(b1,b2)
(b2,b3)
(b2,b3)
(b1,b3,b4)
		
		"""
		
		genes = GeneAnnotation("hg19")
		gene_A1 = Gene("[--A1--]")
		gene_A2 = Gene("[--A2--]")
		gene_A3 = Gene("[--A3--]")
		gene_A4 = Gene("[--A4--]")
		gene_A5 = Gene("[--A5--]")
		gene_A6 = Gene("[--A6--]")
		gene_XX = Gene("X")
		
		genes.add_annotation(gene_A1,"1",10000,13000)
		genes.add_annotation(gene_A2,"1",11500,14500)
		genes.add_annotation(gene_A3,"1",13000,17000)
		genes.add_annotation(gene_A4,"1",15000,18500)
		genes.add_annotation(gene_A5,"1",15000,19000)
		genes.add_annotation(gene_A5,"1",11500,12500)# Add twice
		genes.add_annotation(gene_A6,"1",17000,19000)
		genes.add_annotation(gene_XX,"X",14000,16000)
		
		fusion_1 = Fusion("chr1","chrX",12000,15000,None,None,"+","+","Experiment_1")
		experiment_1 = FusionDetectionExperiment("Experiment_1","RNA")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chrX",14000,15000,None,None,"+","+","Experiment_2")
		experiment_2 = FusionDetectionExperiment("Experiment_2","RNA")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		fusion_3 = Fusion("chr1","chrX",16000,15000,None,None,"+","+","Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3","RNA")
		experiment_3.add_fusion(fusion_3)
		experiment_3.annotate_genes(genes)
		
		fusion_4 = Fusion("chr1","chrX",18000,15000,None,None,"+","+","Experiment_3")
		experiment_4 = FusionDetectionExperiment("Experiment_4","RNA")
		experiment_4.add_fusion(fusion_4)
		experiment_4.annotate_genes(genes)
		
		#experiment_1.show_me()
		#experiment_2.show_me()
		#experiment_3.show_me()
		#experiment_4.show_me()
		
		#(b1,b2) = (A2)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A2--]')
		
		#(b1,b3) = (A5)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b1,b4) = (A2)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b2,b3) = (A3)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A3--]')
		
		#(b2,b4) = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4) = (A4,A5)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 2)
		self.assertTrue( \
			(str(overlap[0][0].annotated_genes_left[0]) == '[--A4--]' and str(overlap[0][0].annotated_genes_left[1]) == '[--A5--]') or \
			(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]' and str(overlap[0][0].annotated_genes_left[1]) == '[--A4--]') \
			)
		
		#(b1,b2),b3 = (A2),(A3,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b3),b2 = (A5),(A2,A3)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b3),b1 = (A3),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2),b4 = (A2),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b4),b2 = (A5),(A2,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b1 = (A4),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b1,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b3,b4),b1 = (A4),(A1,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		
		#(b2,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4),b2 = (A4),(A2,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,"summary",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2,b3,b4) = none
		fh = open('fusions.txt','w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(False,fh,"list",egm=False,strand_specific_matching=True,overlap_based_matching=True)
		self.assertTrue(len(overlap[0]) == 0)
		


def main():
	unittest.main()

if __name__ == '__main__':
	main()
