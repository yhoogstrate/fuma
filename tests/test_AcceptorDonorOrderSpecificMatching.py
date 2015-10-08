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

import unittest,hashlib,os


from fuma.Fusion import AD_DIRECTION_FORWARD
from fuma.Fusion import AD_DIRECTION_REVERSE

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadFusionMap
from fuma.ParseBED import ParseBED
from fuma.OverlapComplex import OverlapComplex
from fuma.Fusion import Fusion
from fuma.FusionDetectionExperiment import FusionDetectionExperiment
from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation
from fuma.CLI import CLI


class TestAcceptorDonorOrderSpecificMatching(unittest.TestCase):
	def test_01(self):
		"""
#1:
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#2:
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]
		"""
		
		args_on = CLI(['--acceptor-donor-order-specific-matching','-f','summary','-s',''])
		args_off = CLI(['-f','summary','-s',''])
		
		gene_A = Gene("A")
		gene_B = Gene("B")
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",10000,20000)
		genes.add_annotation(gene_B,"1",80000,90000)
		
		fusion_1 = Fusion("chr1","chr1",15000,85000,None,None,"+","+","Experiment_1")
		fusion_1.add_location({'left':[fusion_1.get_left_chromosome(), fusion_1.get_left_break_position()], 'right':[fusion_1.get_right_chromosome(), fusion_1.get_right_break_position()], 'id':1, 'dataset':fusion_1.dataset_name })
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chr1",85000,15000,None,None,"+","+","Experiment_2")
		fusion_2.add_location({ 'left':[fusion_2.get_left_chromosome(), fusion_2.get_left_break_position()], 'right':[fusion_2.get_right_chromosome(), fusion_2.get_right_break_position()], 'id':2, 'dataset':fusion_2.dataset_name })
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		self.assertEqual(len(fusion_1.annotated_genes_left) , 1)
		self.assertEqual(len(fusion_2.annotated_genes_left) , 1)
		
		self.assertEqual( fusion_1.acceptor_donor_direction , AD_DIRECTION_FORWARD )
		self.assertEqual( fusion_2.acceptor_donor_direction , AD_DIRECTION_REVERSE )
		
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args_on)
		
		self.assertEqual(len(overlap[0]), 0)
		
		
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args_off)
		
		self.assertEqual(len(overlap[0]), 1)
	
	def test_02(self):
		"""
#1:
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#2:
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]
		"""
		
		args_off = CLI(['-f','summary','-s',''])
		
		gene_A = Gene("A")
		gene_B = Gene("B")
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",10000,20000)
		genes.add_annotation(gene_B,"1",80000,90000)
		
		fusion_1 = Fusion("chr1","chr1",15000,85000,None,None,"+","+","Experiment_1")
		fusion_1.add_location({'left':[fusion_1.get_left_chromosome(), fusion_1.get_left_break_position()], 'right':[fusion_1.get_right_chromosome(), fusion_1.get_right_break_position()], 'id':1, 'dataset':fusion_1.dataset_name })
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		fusion_2 = Fusion("chr1","chr1",85000,15000,None,None,"+","+","Experiment_2")
		fusion_2.add_location({ 'left':[fusion_2.get_left_chromosome(), fusion_2.get_left_break_position()], 'right':[fusion_2.get_right_chromosome(), fusion_2.get_right_break_position()], 'id':2, 'dataset':fusion_2.dataset_name })
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		self.assertEqual(len(fusion_1.annotated_genes_left) , 1)
		self.assertEqual(len(fusion_2.annotated_genes_left) , 1)
		
		self.assertEqual( fusion_1.acceptor_donor_direction , AD_DIRECTION_FORWARD )
		self.assertEqual( fusion_2.acceptor_donor_direction , AD_DIRECTION_REVERSE )
		
		# -> & <- = unknown
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args_off)
		
		self.assertEqual(len(overlap[0]), 1)
		self.assertEqual(overlap[0][0].acceptor_donor_direction, None)
		
		# <- & -> = unknown
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,args_off)
		
		self.assertEqual(len(overlap[0]), 1)
		self.assertEqual(overlap[0][0].acceptor_donor_direction, None)
		
		# -> & -> = ->
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,args_off)
		
		self.assertEqual(len(overlap[0]), 1)
		self.assertEqual(overlap[0][0].acceptor_donor_direction, AD_DIRECTION_FORWARD )
		
		# <- & <- = <-
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args_off)
		
		self.assertEqual(len(overlap[0]), 1)
		self.assertEqual(overlap[0][0].acceptor_donor_direction, AD_DIRECTION_REVERSE )
	
	def test_03(self):
		"""
#AB1:
        -->                       -->
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#AB2:
      <--                       <--
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#AB3:
        -->                     <--
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#AB4:
      <--                         -->
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#BA1:
        -->                       -->
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]

#BA2:
      <--                       <--
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]

#BA3:
        -->                     <--
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]

#BA4:
      <--                         -->
        break1                    break2
        |                         |
[ --- Gene B --- ]        [ --- Gene A --- ]

"""
		gene_A = Gene("A")
		gene_B = Gene("B")
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",10000,20000)
		genes.add_annotation(gene_B,"1",80000,90000)
		
		fusion_AB1 = Fusion("chr1","chr1",15000,85000,None,None,"+","+","Experiment_AB1")
		fusion_AB1.add_location({'left':[fusion_AB1.get_left_chromosome(), fusion_AB1.get_left_break_position()], 'right':[fusion_AB1.get_right_chromosome(), fusion_AB1.get_right_break_position()], 'id':1, 'dataset':fusion_AB1.dataset_name })
		
		fusion_AB2 = Fusion("chr1","chr1",15000,85000,None,None,"-","-","Experiment_AB2")
		fusion_AB2.add_location({'left':[fusion_AB2.get_left_chromosome(), fusion_AB2.get_left_break_position()], 'right':[fusion_AB2.get_right_chromosome(), fusion_AB2.get_right_break_position()], 'id':2, 'dataset':fusion_AB2.dataset_name })
		
		fusion_AB3 = Fusion("chr1","chr1",15000,85000,None,None,"+","-","Experiment_AB3")
		fusion_AB3.add_location({'left':[fusion_AB3.get_left_chromosome(), fusion_AB3.get_left_break_position()], 'right':[fusion_AB3.get_right_chromosome(), fusion_AB3.get_right_break_position()], 'id':3, 'dataset':fusion_AB3.dataset_name })
		
		fusion_AB4 = Fusion("chr1","chr1",15000,85000,None,None,"-","+","Experiment_AB4")
		fusion_AB4.add_location({'left':[fusion_AB4.get_left_chromosome(), fusion_AB4.get_left_break_position()], 'right':[fusion_AB4.get_right_chromosome(), fusion_AB4.get_right_break_position()], 'id':4, 'dataset':fusion_AB4.dataset_name })
		
		
		fusion_BA1 = Fusion("chr1","chr1",85000,15000,None,None,"+","+","Experiment_BA1")
		fusion_BA1.add_location({'left':[fusion_BA1.get_left_chromosome(), fusion_BA1.get_left_break_position()], 'right':[fusion_BA1.get_right_chromosome(), fusion_BA1.get_right_break_position()], 'id':5, 'dataset':fusion_BA1.dataset_name })
		
		fusion_BA2 = Fusion("chr1","chr1",85000,15000,None,None,"-","-","Experiment_BA2")
		fusion_BA2.add_location({'left':[fusion_BA2.get_left_chromosome(), fusion_BA2.get_left_break_position()], 'right':[fusion_BA2.get_right_chromosome(), fusion_BA2.get_right_break_position()], 'id':6, 'dataset':fusion_BA2.dataset_name })
		
		fusion_BA3 = Fusion("chr1","chr1",85000,15000,None,None,"+","-","Experiment_BA3")
		fusion_BA3.add_location({'left':[fusion_BA3.get_left_chromosome(), fusion_BA3.get_left_break_position()], 'right':[fusion_BA3.get_right_chromosome(), fusion_BA3.get_right_break_position()], 'id':7, 'dataset':fusion_BA3.dataset_name })
		
		fusion_BA4 = Fusion("chr1","chr1",85000,15000,None,None,"-","+","Experiment_BA4")
		fusion_BA4.add_location({'left':[fusion_BA4.get_left_chromosome(), fusion_BA4.get_left_break_position()], 'right':[fusion_BA4.get_right_chromosome(), fusion_BA4.get_right_break_position()], 'id':8, 'dataset':fusion_BA4.dataset_name })
		
		experiments = {'AB':[],'BA':[]}
		
		experiments['AB'].append(FusionDetectionExperiment("Experiment_AB1"))
		experiments['AB'][0].add_fusion(fusion_AB1)
		experiments['AB'][0].annotate_genes(genes)
		
		experiments['AB'].append(FusionDetectionExperiment("Experiment_AB2"))
		experiments['AB'][1].add_fusion(fusion_AB2)
		experiments['AB'][1].annotate_genes(genes)
		
		experiments['AB'].append(FusionDetectionExperiment("Experiment_AB3"))
		experiments['AB'][2].add_fusion(fusion_AB3)
		experiments['AB'][2].annotate_genes(genes)
		
		experiments['AB'].append(FusionDetectionExperiment("Experiment_AB4"))
		experiments['AB'][3].add_fusion(fusion_AB4)
		experiments['AB'][3].annotate_genes(genes)
		
		experiments['BA'].append(FusionDetectionExperiment("Experiment_BA1"))
		experiments['BA'][0].add_fusion(fusion_BA1)
		experiments['BA'][0].annotate_genes(genes)
		
		experiments['BA'].append(FusionDetectionExperiment("Experiment_BA2"))
		experiments['BA'][1].add_fusion(fusion_BA2)
		experiments['BA'][1].annotate_genes(genes)
		
		# Swap 3 and 4 - to match "-" , "+" and AB <-> BA
		experiments['BA'].append(FusionDetectionExperiment("Experiment_BA4"))
		experiments['BA'][2].add_fusion(fusion_BA4)
		experiments['BA'][2].annotate_genes(genes)
		
		experiments['BA'].append(FusionDetectionExperiment("Experiment_BA3"))
		experiments['BA'][3].add_fusion(fusion_BA3)
		experiments['BA'][3].annotate_genes(genes)
		
		
		# No strict settings - everything should match with everything
		args = CLI(['-f','summary','-s',''])
		for ad_direction_1 in experiments.keys():
			for breakpoint_strand_1 in range(len(experiments[ad_direction_1])):
				for ad_direction_2 in experiments.keys():
					for breakpoint_strand_2 in range(len(experiments[ad_direction_2])):
						overlapping_complex = OverlapComplex()
						overlapping_complex.add_experiment(experiments[ad_direction_1][breakpoint_strand_1])
						overlapping_complex.add_experiment(experiments[ad_direction_2][breakpoint_strand_2])
						overlap = overlapping_complex.overlay_fusions(True,False,args)
						
						self.assertEqual(len(overlap[0]), 1)
						
						if(ad_direction_1 == ad_direction_2):
							self.assertNotEqual(overlap[0][0].acceptor_donor_direction, None)
						else:
							self.assertEqual(overlap[0][0].acceptor_donor_direction, None)
		
		# No strict settings - everything should match with everything
		args = CLI(['-f','summary','--strand-specific-matching','-s',''])
		for ad_direction_1 in experiments.keys():
			for breakpoint_strand_1 in range(len(experiments[ad_direction_1])):
				for ad_direction_2 in experiments.keys():
					for breakpoint_strand_2 in range(len(experiments[ad_direction_2])):
						overlapping_complex = OverlapComplex()
						overlapping_complex.add_experiment(experiments[ad_direction_1][breakpoint_strand_1])
						overlapping_complex.add_experiment(experiments[ad_direction_2][breakpoint_strand_2])
						overlap = overlapping_complex.overlay_fusions(True,False,args)
						
						if(breakpoint_strand_1 == breakpoint_strand_2):
							self.assertEqual(len(overlap[0]), 1)
						else:
							self.assertEqual(len(overlap[0]), 0)
		
		# No strict settings - everything should match with everything
		args = CLI(['-f','summary','--acceptor-donor-order-specific-matching','-s',''])
		for ad_direction_1 in experiments.keys():
			for breakpoint_strand_1 in range(len(experiments[ad_direction_1])):
				for ad_direction_2 in experiments.keys():
					for breakpoint_strand_2 in range(len(experiments[ad_direction_2])):
						overlapping_complex = OverlapComplex()
						overlapping_complex.add_experiment(experiments[ad_direction_1][breakpoint_strand_1])
						overlapping_complex.add_experiment(experiments[ad_direction_2][breakpoint_strand_2])
						overlap = overlapping_complex.overlay_fusions(True,False,args)
						
						if(ad_direction_1 == ad_direction_2):
							self.assertEqual(len(overlap[0]), 1)
						else:
							self.assertEqual(len(overlap[0]), 0)
		
		# No strict settings - everything should match with everything
		args = CLI(['-f','summary','--strand-specific-matching','--acceptor-donor-order-specific-matching','-s',''])
		for ad_direction_1 in experiments.keys():
			for breakpoint_strand_1 in range(len(experiments[ad_direction_1])):
				for ad_direction_2 in experiments.keys():
					for breakpoint_strand_2 in range(len(experiments[ad_direction_2])):
						overlapping_complex = OverlapComplex()
						overlapping_complex.add_experiment(experiments[ad_direction_1][breakpoint_strand_1])
						overlapping_complex.add_experiment(experiments[ad_direction_2][breakpoint_strand_2])
						overlap = overlapping_complex.overlay_fusions(True,False,args)
						
						if (breakpoint_strand_1 == breakpoint_strand_2) and (ad_direction_1 == ad_direction_2):
							self.assertEqual(len(overlap[0]), 1)
						else:
							self.assertEqual(len(overlap[0]), 0)
	
	def test_04(self):
		"""
#1:
        -->                     <--
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

#2:
      <--                         -->
        break1                    break2
        |                         |
[ --- Gene A --- ]        [ --- Gene B --- ]

Ensure the strand of the merged fusion is not set!
"""
		gene_A = Gene("A")
		gene_B = Gene("B")
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",10000,20000)
		genes.add_annotation(gene_B,"1",80000,90000)
		
		fusion_1 = Fusion("chr1","chr1",15000,85000,None,None,"+","-","Experiment_1")
		fusion_1.add_location({'left':[fusion_1.get_left_chromosome(), fusion_1.get_left_break_position()], 'right':[fusion_1.get_right_chromosome(), fusion_1.get_right_break_position()], 'id':3, 'dataset':fusion_1.dataset_name })
		
		fusion_2 = Fusion("chr1","chr1",15000,85000,None,None,"-","+","Experiment_2")
		fusion_2.add_location({'left':[fusion_2.get_left_chromosome(), fusion_2.get_left_break_position()], 'right':[fusion_2.get_right_chromosome(), fusion_2.get_right_break_position()], 'id':4, 'dataset':fusion_2.dataset_name })
		
		experiment_1 = FusionDetectionExperiment("Experiment_1")
		experiment_1.add_fusion(fusion_1)
		experiment_1.annotate_genes(genes)
		
		experiment_2 = FusionDetectionExperiment("Experiment_2")
		experiment_2.add_fusion(fusion_2)
		experiment_2.annotate_genes(genes)
		
		args = args = CLI(['-f','summary','-s',''])
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args)
		
		self.assertEqual(overlap[0][0].left_strand, None)
		self.assertEqual(overlap[0][0].right_strand, None)
		
		args = args = CLI(['-f','summary','-s',''])
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_1)
		overlapping_complex.add_experiment(experiment_1)
		overlap = overlapping_complex.overlay_fusions(True,False,args)
		
		self.assertNotEqual(overlap[0][0].left_strand, None)
		self.assertNotEqual(overlap[0][0].right_strand, None)
		
		args = args = CLI(['-f','summary','-s',''])
		overlapping_complex = OverlapComplex()
		overlapping_complex.add_experiment(experiment_2)
		overlapping_complex.add_experiment(experiment_2)
		overlap = overlapping_complex.overlay_fusions(True,False,args)
		
		self.assertNotEqual(overlap[0][0].left_strand, None)
		self.assertNotEqual(overlap[0][0].right_strand, None)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
