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

import unittest,hashlib,os,logging,sys
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadFusionMap
from fuma.ParseBED import ParseBED
from fuma.OverlapComplex import OverlapComplex
from fuma.Fusion import Fusion
from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation
from fuma.CLI import CLI

class TestOverlapComplex(unittest.TestCase):
	# The test indicated a minor bug. The unittest aborts on this test
	# case. A solution is given at the following url:
	# <http://stackoverflow.com/questions/4732827/continuing-in-pythons-unittest-when-an-assertion-fails>
	#
	#def setUp(self):
	#	self.verificationErrors = []
	#
	#def tearDown(self):
	#	self.assertEqual([], self.verificationErrors)
	
	def test_01(self):
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_01.txt','w'),args)
		
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
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)#A,B,C
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_1","uid",True)#B,C
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_2","uid",True)#A,B
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_02.txt','w'),args)
		
		#overlapping_complex.export_summary("-")
		
		self.assertEqual(overlapping_complex.matches_total['1.2'],2)
		self.assertEqual(overlapping_complex.matches_total['1.2.3'],0)
	
	def test_04(self):
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_04.txt','w'),args)
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 1)
	
	def test_05(self):
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr3",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_05.txt','w'),args)
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == 1)
		self.assertTrue(overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 0)
	
	def test_06(self):
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,"+","+","Experiment_1","uid",True)
		fusion_2 = Fusion("chrX","chr2",15000,80000,"+","+","Experiment_2","uid",True)
		fusion_3 = Fusion("chrX","chr2",15000,70000,"+","+","Experiment_3","uid",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		
		experiment_1.genes_spanning_left_junction = [True]
		experiment_1.genes_spanning_right_junction = [True]
		experiment_2.genes_spanning_left_junction = [True]
		experiment_2.genes_spanning_right_junction = [True]
		experiment_3.genes_spanning_left_junction = [True]
		experiment_3.genes_spanning_right_junction = [True]
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_06.txt','w'),args)
		
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
		args = CLI(['-m','subset','-f','list','--no-strand-specific-matching','-s',''])
		
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
		
		experiment_1.remove_duplicates(args)
		experiment_2.remove_duplicates(args)
		experiment_3.remove_duplicates(args)
		
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
		
		overlap_1 = overlapping_complex_1.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_1.txt','w'),args)
		overlap_2 = overlapping_complex_2.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_2.txt','w'),args)
		overlap_3 = overlapping_complex_3.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_3.txt','w'),args)
		overlap_4 = overlapping_complex_4.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_4.txt','w'),args)
		overlap_5 = overlapping_complex_5.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_5.txt','w'),args)
		overlap_6 = overlapping_complex_6.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_07_6.txt','w'),args)
		
		self.assertEqual(overlapping_complex_1.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_2.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_3.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_4.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_5.matches_total['1.2.3'],0)
		self.assertEqual(overlapping_complex_6.matches_total['1.2.3'],0)
	
	def test_08(self):
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
		args_a = CLI(['-m','subset','--no-strand-specific-matching','-f','list','-s',''])
		args_b = CLI(['-m','subset',   '--strand-specific-matching','-f','list','-s',''])
		
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
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_08_a.txt','w'),args_a)
		self.assertTrue(len(overlap[0]) == 1)
		
		# Matching, do take strand-specific-matching into account.
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		# FuMa the old way makes a mistake over here
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_08_b.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_b)
		self.assertTrue(len(overlap[0]) == 0)
	
	def test_09(self):
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
		args      = CLI(['-m','overlap','-f','list','--strand-specific-matching','-s',''])
		
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
		
		#(b1,b2) = (A2)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_01.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A2--]')
		
		#(b1,b3) = (A5)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_02.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b1,b4) = (A2)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_03.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b2,b3) = (A3)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_04.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A3--]')
		
		#(b2,b4) = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		# Output incomplete
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_05.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4) = (A4,A5)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_06.output.txt','w'),args)
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
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_07.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b3),b2 = (A5),(A2,A3)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_08.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b3),b1 = (A3),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_09.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2),b4 = (A2),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_10.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b4),b2 = (A5),(A2,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_11.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b1 = (A4),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_12.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_13.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b1,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_14.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b3,b4),b1 = (A4),(A1,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_15.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		
		#(b2,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_16.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_17.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4),b2 = (A4),(A2,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_09_18.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2,b3,b4) = none
		test_filename = 'test_OverlapComplex.TestOverlapComplex.test_09_19.output.txt'
		fh = open(test_filename,'w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(False,fh,args)
		fh.close()
		self.assertTrue(len(overlap[0]) == 0)
		
		md5_input   = hashlib.md5(open(test_filename, 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/'+test_filename, 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove(test_filename)
	
	def test_10(self):
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
		args = CLI(['-m','overlap','-f','list','--strand-specific-matching','-s',''])
		args_summary = CLI(['-m','overlap','-f','summary','--strand-specific-matching','-s',''])
		
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
		
		
		#(b1,b2) = (A2)
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_01.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A2--]')
		
		#(b1,b3) = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_02.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b4) = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_03.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b3) = [--A3--],[--A4--]
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_04.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 2)
		self.assertTrue( \
			(str(overlap[0][0].annotated_genes_left[0]) == '[--A3--]' and str(overlap[0][0].annotated_genes_left[1]) == '[--A4--]') or \
			(str(overlap[0][0].annotated_genes_left[0]) == '[--A4--]' and str(overlap[0][0].annotated_genes_left[1]) == '[--A3--]') \
		)
		
		#(b2,b4) = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_05.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4) = [--A5--]
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_06.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 1)
		self.assertTrue(len(overlap[0][0].annotated_genes_left) == 1)
		self.assertTrue(str(overlap[0][0].annotated_genes_left[0]) == '[--A5--]')
		
		#(b1,b2),b3 = (A2),(A3,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_07.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b3),b2 = (A5),(A2,A3)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_08.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b3),b1 = (A3),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_09.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2),b4 = (A2),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_10.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b4),b2 = (A5),(A2,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_11.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b1 = (A4),(A1,A2,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_12.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_13.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b1,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_14.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4),b1 = (A4),(A1,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_1)
		#open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_15.output.txt','w')
		overlap = overlapping_complex.overlay_fusions(True,False,args_summary)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b2,b3),b4 = (A3),(A4,A4,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_16.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b2,b4),b3 = (A5),(A3,A4)       = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_17.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		#(b3,b4),b2 = (A4),(A2,A3,A5)    = none
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,open('tests/data/test_OverlapComplex.TestOverlapComplex.test_10_18.output.txt','w'),args)
		self.assertTrue(len(overlap[0]) == 0)
		
		
		#(b1,b2,b3,b4) = none
		test_filename = 'test_OverlapComplex.TestOverlapComplex.test_10_19.output.txt'
		fh = open(test_filename,'w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlap = overlapping_complex.overlay_fusions(False,fh,args)
		fh.close()
		self.assertTrue(len(overlap[0]) == 0)
		
		md5_input   = hashlib.md5(open(test_filename, 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/'+test_filename, 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove(test_filename)
	
	def test_11(self):
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
		args = CLI(['-m','overlap','-f','list','--strand-specific-matching','-s',''])
		
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
		
		test_filename = 'test_OverlapComplex.TestOverlapComplex.test_11.output.txt'
		fh = open(test_filename,'w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlapping_complex.add_experiment(experiment_4)
		overlapping_complex.add_experiment(experiment_5)
		overlap = overlapping_complex.overlay_fusions(False,fh,args)
		fh.close()
		self.assertTrue(len(overlap[0]) == 0)
		
		md5_input   = hashlib.md5(open(test_filename, 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/'+test_filename, 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove(test_filename)
	
	def test_12(self):
		"""
		Tests whether the overlap() matching function is implemented correctly 
- This is the two examples given in the github manual -
		"""
		args = CLI(['-m','overlap','-f','summary','-s',''])
		
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
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(False,False,args)
		
		self.assertEqual(len(overlap[0]) , 1)
		self.assertEqual(len(overlap[0][0].annotated_genes_left) , 1)
		self.assertEqual(str(overlap[0][0].annotated_genes_left[0]) , 'BLUE')
		self.assertEqual(str(overlap[0][0].annotated_genes_right[0]) , 'YELLOW')
	
	def test_13(self):
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
		args_subset = CLI(['-m','subset','-f','list','-s',''])
		args_overlap = CLI(['-m','overlap','-f','list','-s',''])
		
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
		
		
		test_filename = 'test_OverlapComplex.TestOverlapComplex.test_13.output-subset.txt'
		fh = open(test_filename,'w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(False,fh,args_subset)
		fh.close()
		self.assertTrue(len(overlap[0]) == 0)
		
		md5_input   = hashlib.md5(open(test_filename, 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/'+test_filename, 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove(test_filename)
		
		
		test_filename = 'test_OverlapComplex.TestOverlapComplex.test_13.output-overlap.txt'
		fh = open(test_filename,'w')
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		overlap = overlapping_complex.overlay_fusions(False,fh,args_overlap)
		fh.close()
		self.assertTrue(len(overlap[0]) == 1)
		
		md5_input   = hashlib.md5(open(test_filename, 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/'+test_filename, 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove(test_filename)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
