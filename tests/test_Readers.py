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

import sys,unittest,logging
logging.basicConfig(level=logging.INFO,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Fusion import STRAND_FORWARD
from fuma.Fusion import STRAND_REVERSE

from fuma.Fusion import AD_DIRECTION_FORWARD
from fuma.Fusion import AD_DIRECTION_REVERSE


from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadChimeraPrettyPrint
from fuma.Readers import ReadRNASTARFusionFinal
from fuma.Readers import ReadSOAPFuseGenes
from fuma.Readers import ReadSOAPFuseTranscripts
from fuma.Readers import ReadEricScriptResultsTotal
from fuma.Readers import ReadJaffaResults


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
		
		self.assertEqual(fusions[0].left_break_position , 8407 )
		self.assertEqual(fusions[0].right_break_position , 14006 )
		self.assertEqual(fusions[0].left_strand , STRAND_REVERSE )
		self.assertEqual(fusions[0].right_strand , STRAND_FORWARD )
		
		self.assertEqual(fusions[1].left_break_position , 5320 )
		self.assertEqual(fusions[1].right_break_position , 11706 )
		self.assertEqual(fusions[1].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[1].right_strand , STRAND_REVERSE )
	
	def test_02(self):
		""" Tests whether files of input format Chimera prettyPrint can
		be parsed"""
		
		fusions = ReadChimeraPrettyPrint("tests/data/test_Readers.TestReadChimeraPrettyPrint.test_02.txt","test")
		
		self.assertEqual(fusions[0].acceptor_donor_direction , AD_DIRECTION_FORWARD )
		self.assertEqual(fusions[1].acceptor_donor_direction , AD_DIRECTION_REVERSE )
		self.assertEqual(fusions[2].acceptor_donor_direction , AD_DIRECTION_REVERSE )
		self.assertEqual(fusions[3].acceptor_donor_direction , AD_DIRECTION_FORWARD )


class TestReadRNASTARFusionFinal(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from Star Fusion can
		be parsed"""
		
		fusions = ReadRNASTARFusionFinal("tests/data/test_Readers.TestReadRNASTARFusionFinal.test_01.candidates.final","test")
		
		self.assertEqual(len(fusions) , 3)
		
		# 'chr11' < 'chr5' -> swap left and right
		self.assertEqual(fusions[0].left_break_position , 77580768 )
		self.assertEqual(fusions[0].right_break_position , 1085347 )
		self.assertEqual(fusions[0].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[0].right_strand , STRAND_REVERSE )
		
		self.assertEqual(fusions[1].left_break_position , 156294763 )
		self.assertEqual(fusions[1].right_break_position , 156374393 )
		self.assertEqual(fusions[1].left_strand , STRAND_REVERSE )
		self.assertEqual(fusions[1].right_strand , STRAND_REVERSE )
		
		self.assertEqual(fusions[2].left_break_position , 29446010 )
		self.assertEqual(fusions[2].right_break_position , 43479658 )
		self.assertEqual(fusions[2].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[2].right_strand , STRAND_REVERSE )


class TestCompleteGenomics(unittest.TestCase):
	def test_01(self):
		#self.failUnless(1 == 1)
		#self.failIf(1 != 1)
		pass


class TestReadSOAPFuseGenes(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from SOAPFusion can
		be parsed"""
		
		fusions = ReadSOAPFuseGenes("tests/data/test_Readers.TestReadSOAPFuseGenes.test_01.txt","test")
		
		self.assertEqual(len(fusions) , 3)
		
		self.assertEqual(fusions[0].get_left_chromosome(True) , 'chr11')
		self.assertEqual(fusions[0].get_right_chromosome(True) , 'chr9')
		self.assertEqual(fusions[0].left_break_position , 65267181 )
		self.assertEqual(fusions[0].right_break_position , 32542249 )
		self.assertEqual(fusions[0].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[0].right_strand , STRAND_REVERSE )
		
		self.assertEqual(fusions[1].get_left_chromosome(True) , 'chr11')
		self.assertEqual(fusions[1].get_right_chromosome(True) , 'chr19')
		self.assertEqual(fusions[1].left_break_position , 125479335 )
		self.assertEqual(fusions[1].right_break_position , 55899383 )
		self.assertEqual(fusions[1].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[1].right_strand , STRAND_FORWARD )
		
		self.assertEqual(fusions[2].get_left_chromosome(True) , 'chr11')
		self.assertEqual(fusions[2].get_right_chromosome(True) , 'chr20')
		self.assertEqual(fusions[2].left_break_position , 65116155 )
		self.assertEqual(fusions[2].right_break_position , 33114078 )
		self.assertEqual(fusions[2].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[2].right_strand , STRAND_FORWARD )


class TestReadSOAPFuseTranscripts(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from SOAPFusion can
		be parsed"""
		
		fusions = ReadSOAPFuseTranscripts("tests/data/test_Readers.TestReadSOAPFuseTranscripts.test_01.txt","test")
		
		self.assertEqual(len(fusions) , 3)
		
		self.assertEqual(fusions[0].get_left_chromosome(True) , 'chr11')
		self.assertEqual(fusions[0].get_right_chromosome(True) , 'chr20')
		self.assertEqual(fusions[0].left_break_position , 65116155 )
		self.assertEqual(fusions[0].right_break_position , 33114078 )
		self.assertEqual(fusions[0].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[0].right_strand , STRAND_FORWARD )
		
		self.assertEqual(fusions[1].get_left_chromosome(True) , 'chr11')
		self.assertEqual(fusions[1].get_right_chromosome(True) , 'chr20')
		self.assertEqual(fusions[1].left_break_position , 65116155 )
		self.assertEqual(fusions[1].right_break_position , 33114078 )
		self.assertEqual(fusions[1].left_strand , STRAND_FORWARD )
		self.assertEqual(fusions[1].right_strand , STRAND_FORWARD )
		
		self.assertEqual(fusions[2].get_left_chromosome(True) , 'chr7')
		self.assertEqual(fusions[2].get_right_chromosome(True) , 'chr7')
		self.assertEqual(fusions[2].left_break_position , 99521226 )
		self.assertEqual(fusions[2].right_break_position , 99569369 )
		self.assertEqual(fusions[2].left_strand , STRAND_REVERSE )
		self.assertEqual(fusions[2].right_strand , STRAND_REVERSE )


class TestReadEricScriptResultsTotal(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from SOAPFusion can
		be parsed"""
		
		fusions = ReadEricScriptResultsTotal("tests/data/test_Readers.TestReadEricScriptResultsTotal.test_01.txt","test")
		
		self.assertEqual(len(fusions) , 2)
		
		self.assertEqual(fusions[0].get_left_chromosome(True) , 'chr1')
		self.assertEqual(fusions[0].get_right_chromosome(True) , 'chrX')
		self.assertEqual(fusions[0].left_break_position , 94976043 )
		self.assertEqual(fusions[0].right_break_position , 123911702 )
		self.assertEqual(fusions[0].left_strand , STRAND_REVERSE )
		self.assertEqual(fusions[0].right_strand , STRAND_FORWARD )
		
		self.assertEqual(fusions[1].get_left_chromosome(True) , 'chr3')
		self.assertEqual(fusions[1].get_right_chromosome(True) , 'chrX')
		self.assertEqual(fusions[1].left_break_position , 156540631 )
		self.assertEqual(fusions[1].right_break_position , 123911171 )
		self.assertEqual(fusions[1].left_strand , STRAND_REVERSE )
		self.assertEqual(fusions[1].right_strand , STRAND_FORWARD )


class TestReadJaffaResults(unittest.TestCase):
	def test_01(self):
		""" Tests whether files of input format from SOAPFusion can
		be parsed
		
		Fusion (from dataset 'test'): chr20:46365686(?)<->chr20:47538547(?)
		Fusion (from dataset 'test'): chr17:59445688(?)<->chr20:49411710(?)
		Fusion (from dataset 'test'): chr17:37793484(?)<->chr20:53259997(?)
		Fusion (from dataset 'test'): chr17:57917129(?)<->chr17:57992064(?)
		"""
		
		fusions = ReadJaffaResults("tests/data/test_Readers.TestReadJaffaResults.test_01.txt","test")
		
		self.assertEqual(len(fusions) , 4)
		
		self.assertEqual(fusions[0].get_left_chromosome(True) , 'chr20')
		self.assertEqual(fusions[0].get_right_chromosome(True) , 'chr20')
		self.assertEqual(fusions[0].left_break_position , 46365686)
		self.assertEqual(fusions[0].right_break_position , 47538547)
		self.assertEqual(fusions[0].left_strand , None)
		self.assertEqual(fusions[0].right_strand , None)
		self.assertEqual(fusions[0].acceptor_donor_direction , None)
		
		#dataset 'test'): chr17:59445688(?)<-chr20:49411710(?)
		self.assertEqual(fusions[1].get_left_chromosome(True) , 'chr17')
		self.assertEqual(fusions[1].get_right_chromosome(True) , 'chr20')
		self.assertEqual(fusions[1].left_break_position , 59445688)
		self.assertEqual(fusions[1].right_break_position , 49411710)
		self.assertEqual(fusions[1].left_strand , None)
		self.assertEqual(fusions[1].right_strand , None)
		self.assertEqual(fusions[0].acceptor_donor_direction , None)
		
		# @todo
		# comparing 2x test read jaffa results should give a exception:
		# raise Exception("A fusion gene without an annotated acceptor-donor direction was used for acceptor-donor-order-specific-matching.\n\n"+fusion_1.__str__()+"\n"+fusion_2.__str__())


def main():
	unittest.main()

if __name__ == '__main__':
	main()
