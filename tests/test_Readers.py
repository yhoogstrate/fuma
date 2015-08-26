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

import sys,unittest

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadChimeraPrettyPrint
from fuma.Readers import ReadRNASTARFusionFinal


class TestReadChimeraScanAbsoluteBEDPE(unittest.TestCase):
	def test_01(self):
		fusions = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Readers.TestReadChimeraScanAbsoluteBEDPE.test_01.bedpe","RNA")
		i = 0
		
		for fusion in fusions:
			i += 1
		
		self.assertEqual(i, 690)
	
	def test_02(self):
		"""
		In this test we have the same breakpoint twice, but swapped
		the 5' and 3'. After a duplicaiton removal only one fusion
		should be kept.
		"""
		fusions = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Readers.TestReadChimeraScanAbsoluteBEDPE.test_02.bedbe","RNA")
		
		j = 0
		
		for fusion in fusions:
			left_break = fusion.get_left_position(True)
			right_break = fusion.get_right_position(True)
			
			self.failUnless(left_break[0] == 'chr4')
			self.assertEqual(left_break[1] , 77000000)
			
			self.failUnless(right_break[0] == 'chr7')
			self.assertEqual(right_break[1] , 20000000)
			
			j += 1
		
		self.failUnless(j == 4)

class TestReadChimeraPrettyPrint(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format Chimera prettyPrint can
		be parsed"""
		
		fusions = ReadChimeraPrettyPrint("tests/data/test_Readers.TestReadChimeraPrettyPrint.test_01.txt","test")
		
		self.assertEqual(len(fusions) , 2)
		
		self.assertEqual( fusions[0].left_break_position , 8407 )
		self.assertEqual( fusions[0].right_break_position , 14006 )
		self.assertEqual( fusions[0].left_strand , STRAND_REVERSE )
		self.assertEqual( fusions[0].right_strand , STRAND_FORWARD )
		
		self.assertEqual( fusions[1].left_break_position , 5320 )
		self.assertEqual( fusions[1].right_break_position , 11706 )
		self.assertEqual( fusions[1].left_strand , STRAND_FORWARD )
		self.assertEqual( fusions[1].right_strand , STRAND_REVERSE )
		
		fusions.show_me()

class TestReadRNASTARFusionFinal(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from Star Fusion can
		be parsed"""
		
		fusions = ReadRNASTARFusionFinal("tests/data/test_Readers.TestReadRNASTARFusionFinal.test_01.candidates.final","test")
		
		self.assertEqual(len(fusions) , 3)
		
		# 'chr11' < 'chr5' -> swap left and right
		self.assertEqual( fusions[0].left_break_position , 77580768 )
		self.assertEqual( fusions[0].right_break_position , 1085347 )
		self.assertEqual( fusions[0].left_strand , STRAND_FORWARD )
		self.assertEqual( fusions[0].right_strand , STRAND_REVERSE )
		
		self.assertEqual( fusions[1].left_break_position , 156294763 )
		self.assertEqual( fusions[1].right_break_position , 156374393 )
		self.assertEqual( fusions[1].left_strand , STRAND_REVERSE )
		self.assertEqual( fusions[1].right_strand , STRAND_REVERSE )
		
		self.assertEqual( fusions[2].left_break_position , 29446010 )
		self.assertEqual( fusions[2].right_break_position , 43479658 )
		self.assertEqual( fusions[2].left_strand , STRAND_FORWARD )
		self.assertEqual( fusions[2].right_strand , STRAND_REVERSE )
		
		fusions.show_me()



class TestCompleteGenomics(unittest.TestCase):
	def test_01(self):
		#self.failUnless(1 == 1)
		#self.failIf(1 != 1)
		pass

def main():
	unittest.main()

if __name__ == '__main__':
	main()
