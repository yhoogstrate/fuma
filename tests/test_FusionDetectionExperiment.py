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

import unittest,logging,sys
logging.basicConfig(level=logging.INFO,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.ParseBED import ParseBED
from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Fusion import Fusion
from fuma.Gene import Gene
from fuma.CLI import CLI

class TestFusionDetectionExperiment(unittest.TestCase):
	def test_01(self):
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		experiment = ReadChimeraScanAbsoluteBEDPE("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_01.bedpe","TestExperiment")
		genes =                          ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_01.bed","hg18", 200000)
		
		length_before_duplication_removal = len(experiment)
		
		experiment.annotate_genes(genes)
		experiment.remove_duplicates(args)
		
		length_after_duplication_removal = len(experiment)
		
		self.assertTrue(length_before_duplication_removal > length_after_duplication_removal)
	
	def test_02(self):
		experiment = ReadChimeraScanAbsoluteBEDPE("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_02.bedpe","TestExperiment")
		genes =                          ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_02.bed","hg18", 200000)
		
		self.assertEqual(len(experiment), 1)
		
		experiment.annotate_genes(genes)
		
		for fusion in experiment:
			left_genes = fusion.get_annotated_genes_left()
			right_genes = fusion.get_annotated_genes_right()
			
			self.assertEqual(len(left_genes), 8)
			self.assertEqual(len(right_genes), 8)
			
			# Ensure all annotated gene names do NOT contain substring '_invalid'
			self.assertEqual(min([str(gene_name).find("_invalid") for gene_name in left_genes]) , -1)
			self.assertEqual(min([str(gene_name).find("_invalid") for gene_name in right_genes]) , -1)
		
	def test_03(self):
		"""
		Check the duplication removal - simple test; 2 identical fusions
		"""
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_03.bed","hg18", 200000)
		experiment.annotate_genes(genes)
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 1)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 1)
			self.assertEqual(len(fusion.annotated_genes_right), 1)
	
	def test_04(self):
		"""
		Check the duplication removal - simple test; 2 identical fusions but checking presevation of the gene names from different annotations
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		gene_1_hg18 = Gene("gene_1", False)
		gene_1_hg19 = Gene("gene_1", False)
		
		gene_2_hg18 = Gene("gene_2", False)
		gene_2_hg19 = Gene("gene_2", False)
		
		fusion_hg18 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_hg19 = Fusion("chr1","chr2",15500,20500,None,None,"+","+","Experiment","",True)
		
		fusion_hg18.annotate_genes_left([gene_1_hg18])
		fusion_hg19.annotate_genes_left([gene_1_hg19])
		
		fusion_hg18.annotate_genes_right([gene_2_hg18])
		fusion_hg19.annotate_genes_right([gene_2_hg19])
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.genes_spanning_left_junction = [True]
		experiment.genes_spanning_right_junction = [True]
		
		experiment.add_fusion(fusion_hg18)
		experiment.add_fusion(fusion_hg19)
		
		self.assertEqual(len(experiment), 2)
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 1)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 2)
			#self.assertEqual(len(fusion.annotated_genes_right), 2)
	
	def test_05(self):
		"""
		Check the duplication removal - 2 fusions; one has two
		overlapping left genes, the other only one, and for right is
		vice versa. The 1 gene is always a subset of the others 2 genes
		and must therefore be treated as identical annotations
		"""
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20050,None,None,"+","+","Experiment","",True)#(1A):(2A,2B)
		fusion_2 = Fusion("chr1","chr2",15050,20000,None,None,"+","+","Experiment","",True)#(1A,1B):(2A)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_05.bed","hg18", 200000)
		experiment.annotate_genes(genes)
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 1)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 1)		# subset = (1A)
			self.assertEqual(len(fusion.annotated_genes_right), 1)		# subset = (2A)
	
	def test_06(self):
		"""
		Check the duplication removal - 2 fusions with only the left
		gene identical
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,30000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_06.bed","hg18", 200000)
		experiment.annotate_genes(genes)
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 2)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 1)
			self.assertEqual(len(fusion.annotated_genes_right), 1)
	
	def test_07(self):
		"""
		Check the duplication removal - 2 fusions with no identical
		genes
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",25000,30000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_07.bed","hg18", 200000)
		
		experiment.annotate_genes(genes)
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 2)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 1)
			self.assertEqual(len(fusion.annotated_genes_right), 1)
	
	def test_08(self):
		"""
		Check the duplication removal - 2 fusions where one is missing
		annotations -> one should be lost because it isn't gene
		spanning
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,30000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_08.bed","hg18", 200000)
		experiment.annotate_genes(genes)
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 1)
		
		for fusion in experiment:
			self.assertEqual(len(fusion.annotated_genes_left), 1)
			self.assertEqual(len(fusion.annotated_genes_right), 1)
	
	def test_09(self):
		"""
		Check the duplication removal - 2 fusions where both are missing
		annotations -> bothshould be lost because it isn't gene
		spanning
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,20000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_09.bed","hg18", 200000)
		experiment.annotate_genes(genes)
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 0)
	
	def test_10(self):
		"""
		May not have overlap:
		{1} = a,b,c
		{2} = a,d
		
		Break {1}:
		                    |
		Break {2}:          :
		       |            :
		a: [   :            :        ]
		b:     :  [         :        ]
		c:     :     [      :        ]
		d: [   : ]          :
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,30000,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,25000,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_10.bed","hg18", 200000)
		
		experiment.annotate_genes(genes)
		
		for fusion in experiment:
			fusion.show_me()
			if(fusion.get_right_break_position() == 30000):
				self.assertEqual(len(fusion.annotated_genes_right), 3)
				self.assertTrue("NM_00002A" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002B" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002C" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
			if(fusion.get_right_break_position() == 25000):
				self.assertEqual(len(fusion.annotated_genes_right), 2)
				self.assertTrue("NM_00002A" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002D" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 2)
	
	def test_11(self):
		"""
		May not have overlap:
		{1} = a,b,c
		{2} = a,b,d
		
		Break {1}:
		                    |
		Break {2}:          :
		       |            :
		a: [   :            :        ]
		b:   [ |            :  ]
		c:     :     [      :        ]
		d: [   : ]     
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr2",15000,27500,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr1","chr2",15000,25500,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_2)
		
		self.assertEqual(len(experiment), 2)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_11.bed","hg18", 200000)
		
		experiment.annotate_genes(genes)
		
		for fusion in experiment:
			if(fusion.get_right_break_position() == 30000):
				self.assertEqual(len(fusion.annotated_genes_right), 3)
				self.assertTrue("NM_00002A" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002B" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002C" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
			if(fusion.get_right_break_position() == 25000):
				self.assertEqual(len(fusion.annotated_genes_right), 3)
				self.assertTrue("NM_00002A" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002B" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
				self.assertTrue("NM_00002D" in [str(gene_name) for gene_name in fusion.annotated_genes_right])
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 2)
	
	def test_12(self):
		"""
		{a} = 1 , 4 , 5
		{b} = 2 , 6
		{c} = 7 , 8
		{d} = 3
		
		   | f1 | f2 | f3 | f4 | f5 | f6 | f7 | f8 |
		---+----+----+----+----+----+----+----+----+
		f1 | *  |    |    | a  | a  |    |    |    |
		---+----+----+----+----+----+----+----+----+
		f2 |    | *  |    |    |    | b  |    |    |
		---+----+----+----+----+----+----+----+----+
		f3 |    |    | d  |    |    |    |    |    |
		---+----+----+----+----+----+----+----+----+
		f4 | a  |    |    | *  | a  |    |    |    |
		---+----+----+----+----+----+----+----+----+
		f5 | a  |    |    | a  | *  |    |    |    |
		---+----+----+----+----+----+----+----+----+
		f6 |    | b  |    |    |    | *  |    |    |
		---+----+----+----+----+----+----+----+----+
		f7 |    |    |    |    |    |    | *  | d  |
		---+----+----+----+----+----+----+----+----+
		f8 |    |    |    |    |    |    | c  | *  |
		---+----+----+----+----+----+----+----+----+
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1 = Fusion("chr1","chr1",15010,80040,None,None,"+","+","Experiment","",True)
		fusion_2 = Fusion("chr2","chr2",15030,80030,None,None,"+","+","Experiment","",True)
		fusion_3 = Fusion("chr4","chr4",15050,80070,None,None,"+","+","Experiment","",True)
		fusion_4 = Fusion("chr1","chr1",15060,80010,None,None,"+","+","Experiment","",True)
		fusion_5 = Fusion("chr1","chr1",15020,80050,None,None,"+","+","Experiment","",True)
		fusion_6 = Fusion("chr2","chr2",15080,80080,None,None,"+","+","Experiment","",True)
		fusion_7 = Fusion("chr3","chr3",15040,80020,None,None,"+","+","Experiment","",True)
		fusion_8 = Fusion("chr3","chr3",15070,80060,None,None,"+","+","Experiment","",True)
		
		experiment = FusionDetectionExperiment("Experiment_1")
		
		experiment.add_fusion(fusion_2)
		experiment.add_fusion(fusion_8)
		experiment.add_fusion(fusion_6)
		experiment.add_fusion(fusion_4)
		experiment.add_fusion(fusion_5)
		experiment.add_fusion(fusion_3)
		experiment.add_fusion(fusion_1)
		experiment.add_fusion(fusion_7)
		
		self.assertEqual(len(experiment), 8)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_12.bed","hg18", 200000)
		
		experiment.annotate_genes(genes)
		
		experiment.remove_duplicates(args)
		
		self.assertEqual(len(experiment), 4)
	
	def test_13(self):
		"""
		Crazy stuff:
		
		{a} = 1 , 2
		{b} = 1 , 3
		
		   | f1 | f2 | f3 |
		---+----+----+----+
		f1 | *  | 1  | 2  |
		---+----+----+----+
		f2 | 1  | *  |    |
		---+----+----+----+
		f3 | 2  |    | *  |
		---+----+----+----+
		
		1 = [a,b,c]
		2 = [a,c]
		3 = [a,b]
		
		>> Try in all possible orders
		
				Break {1}:
		                  |
		Break {2}:        :
		           |      :
		Break {3}: :      :
		           :      :    | 
		a: [       :      :    :     ]
		b:         :  [   :    :     ]
		c: [       :      :  ] :
		
		"""
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		fusion_1_exp_1 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_1","",True)
		fusion_2_exp_1 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_1","",True)
		fusion_3_exp_1 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_1","",True)
		
		fusion_1_exp_2 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_2","",True)
		fusion_2_exp_2 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_2","",True)
		fusion_3_exp_2 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_2","",True)
		
		fusion_1_exp_3 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_3","",True)
		fusion_2_exp_3 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_3","",True)
		fusion_3_exp_3 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_3","",True)
		
		fusion_1_exp_4 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_4","",True)
		fusion_2_exp_4 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_4","",True)
		fusion_3_exp_4 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_4","",True)
		
		fusion_1_exp_5 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_5","",True)
		fusion_2_exp_5 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_5","",True)
		fusion_3_exp_5 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_5","",True)
		
		fusion_1_exp_6 = Fusion("chr1","chr2",15000,70000,None,None,"+","+","Experiment_6","",True)
		fusion_2_exp_6 = Fusion("chr1","chr2",15000,80000,None,None,"+","+","Experiment_6","",True)
		fusion_3_exp_6 = Fusion("chr1","chr2",15000,60000,None,None,"+","+","Experiment_6","",True)
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_3 = FusionDetectionExperiment("Experiment_3")
		experiment_4 = FusionDetectionExperiment("Experiment_4")
		experiment_5 = FusionDetectionExperiment("Experiment_5")
		experiment_6 = FusionDetectionExperiment("Experiment_6")
		
		experiment_1.add_fusion(fusion_1_exp_1)
		experiment_1.add_fusion(fusion_2_exp_1)
		experiment_1.add_fusion(fusion_3_exp_1)
		
		experiment_2.add_fusion(fusion_1_exp_2)
		experiment_2.add_fusion(fusion_3_exp_2)
		experiment_2.add_fusion(fusion_2_exp_2)
		
		experiment_3.add_fusion(fusion_2_exp_3)
		experiment_3.add_fusion(fusion_1_exp_3)
		experiment_3.add_fusion(fusion_3_exp_3)
		
		experiment_4.add_fusion(fusion_2_exp_4)
		experiment_4.add_fusion(fusion_3_exp_4)
		experiment_4.add_fusion(fusion_1_exp_4)
		
		experiment_5.add_fusion(fusion_3_exp_5)
		experiment_5.add_fusion(fusion_1_exp_5)
		experiment_5.add_fusion(fusion_2_exp_5)
		
		experiment_6.add_fusion(fusion_3_exp_6)
		experiment_6.add_fusion(fusion_2_exp_6)
		experiment_6.add_fusion(fusion_1_exp_6)
		
		self.assertEqual(len(experiment_1), 3)
		self.assertEqual(len(experiment_2), 3)
		self.assertEqual(len(experiment_3), 3)
		self.assertEqual(len(experiment_4), 3)
		self.assertEqual(len(experiment_5), 3)
		self.assertEqual(len(experiment_6), 3)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_13.bed","hg18", 200000)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		experiment_4.annotate_genes(genes)
		experiment_5.annotate_genes(genes)
		experiment_6.annotate_genes(genes)
		
		experiment_1.remove_duplicates(args)
		experiment_2.remove_duplicates(args)
		experiment_3.remove_duplicates(args)
		experiment_4.remove_duplicates(args)
		experiment_5.remove_duplicates(args)
		experiment_6.remove_duplicates(args)
		
		# Removing duplicates:
		# 
		# order 1:
		# -> [A,B]   <-> [A,C]   = [A,B]  & [A,C] & [A,B,C]
		# -> [A,B]   <-> [A,B,C] = [A,B*] & [A,C]
		# -> [A,B*]  <-> [A,C]   = [A,B*] & [A,C]
		
		# order 2:
		# -> [A,B]   <-> [A,B,C] = [A,B*] & [A,C]
		# -> [A,B*]  <-> [A,C]   = [A,B*] & [A,C]
		
		# order 3:
		# -> [A,C]   <-> [A,B]   = [A,C]  & [A,B] & [A,B,C]
		# -> [A,C]   <-> [A,B,C] = [A,C*] & [A,B]
		# -> [A,C*]  <-> [A,B]   = [A,C*] & [A,B]
		
		# order 4:
		# -> [A,C]   <-> [A,B,C] = [A,C*] & [A,B]
		# -> [A,C*]  <-> [A,B]   = [A,C*] & [A,B]
		
		# order 5:
		# -> [A,B,C] <-> [A,B]   = [A,B*] & [A,C]
		# -> [A,B*]  <-> [A,C]   = [A,B*] & [A,C]
		
		# order 6:
		# -> [A,B,C] <-> [A,C]   = [A,C*] & [A,B]
		# -> [A,C*]  <-> [A,B]   = [A,C*] & [A,B]
		
		self.assertEqual(len(experiment_1), 2)
		self.assertEqual(len(experiment_2), 2)
		self.assertEqual(len(experiment_3), 2)
		self.assertEqual(len(experiment_4), 2)
		self.assertEqual(len(experiment_5), 2)
		self.assertEqual(len(experiment_6), 2)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
