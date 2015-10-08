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

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE

from fuma.Fusion import AD_DIRECTION_FORWARD
from fuma.Fusion import AD_DIRECTION_REVERSE


from fuma.ParseBED import ParseBED
from fuma.Fusion import Fusion
from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation
from fuma.GeneAnnotation import GeneAnnotation

class TestFusion(unittest.TestCase):
	def test_01(self):
		"""
		gene_A = Gene("A")
		gene_B = Gene("B")
		gene_C = Gene("C")
		gene_X = Gene("X")
		
		genes = GeneAnnotation("hg19")
		genes.add_annotation(gene_A,"1",12000,16000)
		genes.add_annotation(gene_B,"1",13000,17000)
		genes.add_annotation(gene_C,"1",14000,18000)
		genes.add_annotation(gene_X,"X",10000,20000)
		"""
		
		fusion_1 = Fusion("chr1","chrX",15000,15000,None,None,"-","+","Experiment_1")
		fusion_1.add_location({'left':[fusion_1.get_left_chromosome(), fusion_1.get_left_break_position()], 'right':[fusion_1.get_right_chromosome(), fusion_1.get_right_break_position()], 'id':1, 'dataset':fusion_1.dataset_name })
		
		self.assertEqual( fusion_1.left_break_position , 15000 )
		self.assertEqual( fusion_1.right_break_position , 15000 )
		self.assertEqual( fusion_1.left_strand , STRAND_REVERSE )
		self.assertEqual( fusion_1.right_strand , STRAND_FORWARD )

def main():
	unittest.main()

if __name__ == '__main__':
	main()
