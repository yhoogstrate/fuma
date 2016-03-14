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
		fusion_1 = Fusion("chr1","chrX",15000,15000,None,None,"-","+","Experiment_1","1",True)
		
		self.assertEqual( fusion_1.left_break_position , 15000 )
		self.assertEqual( fusion_1.right_break_position , 15000 )
		self.assertEqual( fusion_1.left_strand , STRAND_REVERSE )
		self.assertEqual( fusion_1.right_strand , STRAND_FORWARD )

def main():
	unittest.main()

if __name__ == '__main__':
	main()
