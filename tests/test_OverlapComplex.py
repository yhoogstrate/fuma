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
from fuma.ParseBED import ParseBED
from fuma.OverlapComplex import OverlapComplex
from fuma.Fusion import Fusion
from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Gene import Gene

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
		
		self.assertTrue(len(experiment_1) == len(experiment_2) == len(experiment_3) == 690)
		
		genes = ParseBED("tests/data/test_FusionDetectionExperiment.TestFusionDetectionExperiment.test_01.bed","hg18")
		
		self.assertEqual(len(genes), 47790)
		
		experiment_1.annotate_genes(genes)
		experiment_2.annotate_genes(genes)
		experiment_3.annotate_genes(genes)
		
		experiment_1.remove_duplicates("by-gene-names")
		experiment_2.remove_duplicates("by-gene-names")
		experiment_3.remove_duplicates("by-gene-names")
		
		self.assertTrue(len(experiment_1) == len(experiment_2) == len(experiment_3) <= 690)
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_3)
		
		overlap = overlapping_complex.overlay_fusions()
		
		#overlapping_complex.export_summary("-")
		
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
		f1: [X] -> [A,B,C]
		
		n overlap = 1 (overlap is measured as the merged fusions - only
		one merged fusion can account for all of them)
		"""
		
		fusion_1 = Fusion("chrX","chr2",15000,60000,None,None,"+","+","Experiment_1")
		fusion_2 = Fusion("chrX","chr2",15000,80000,None,None,"+","+","Experiment_1")
		fusion_3 = Fusion("chrX","chr2",15000,70000,None,None,"+","+","Experiment_2")
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
		
		overlap = overlapping_complex.overlay_fusions()
		
		#overlapping_complex.export_summary("-")
		
		self.assertEqual(overlapping_complex.matches_total['1.2'],1)
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
		
		overlap = overlapping_complex.overlay_fusions()
		
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
		
		overlap = overlapping_complex.overlay_fusions()
		
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
		
		overlap = overlapping_complex.overlay_fusions()
		
		overlapping_complex.export_summary("-")
		
		self.assertTrue(overlapping_complex.matches_total['1'] == overlapping_complex.matches_total['2'] == overlapping_complex.matches_total['3'] == overlapping_complex.matches_total['1.2'] == 1)
		self.assertTrue(overlapping_complex.matches_total['1.3'] == overlapping_complex.matches_total['2.3'] == overlapping_complex.matches_total['1.2.3'] == 0)
	
	def test_10(self):
		"""
		Exp1:
		f1: [X] -> [A,B]
		
		Exp2:
		f2: [X] -> [B,C]
		
		Exp3:
		f1: [X] -> [A,B,C]
		
		Expected Exp 1+2+3:
		f1: [X] -> [A,B,C]
		
		n overlap = 1
		
		--------------------------------
		
		@bug depending on the order, the results may vary
		
		@todo this is a very complex use cases and would require a non-
		oob based comparison method. Since the solution is expected to
		be very rare, fixing this issue is not the primary focus
		
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
		
		genes = ParseBED("tests/data/test_OverlapComplex.TestOverlapComplex.test_10.bed","hg18")
		
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
		
		overlap_1 = overlapping_complex_1.overlay_fusions()
		overlap_2 = overlapping_complex_2.overlay_fusions()
		overlap_3 = overlapping_complex_3.overlay_fusions()
		overlap_4 = overlapping_complex_4.overlay_fusions()
		overlap_5 = overlapping_complex_5.overlay_fusions()
		overlap_6 = overlapping_complex_6.overlay_fusions()
		
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_1.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_2.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_3.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_4.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_5.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		#overlapping_complex_6.export_summary("-")
		#print "\n-------------------------------------------------------\n"
		
		try:
			self.assertEqual(overlapping_complex_1.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))
		
		try:
			self.assertEqual(overlapping_complex_2.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))
		
		try:
			self.assertEqual(overlapping_complex_3.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))
		
		try:
			self.assertEqual(overlapping_complex_4.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))
		
		try:
			self.assertEqual(overlapping_complex_5.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))
		
		try:
			self.assertEqual(overlapping_complex_6.matches_total['1.2.3'],1)
		except AssertionError, e: print ("ERROR: "+str(e))

def main():
	unittest.main()

if __name__ == '__main__':
	main()
