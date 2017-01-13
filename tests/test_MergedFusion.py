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



import fuma
from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.ParseBED import ParseBED
from fuma.CLI import CLI
from fuma.MergedFusion import MergedFusion


class TestMergedFusion(unittest.TestCase):
	def test_01(self):
		args = CLI(['-m','egm','--no-strand-specific-matching','-s','','-o','test_ComparisonTriangle.test_02.output.txt'])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		#F[a] + F[b] = MF(A, B)
		mf_a = MergedFusion()
		for fusion in experiment_a:
			mf_a.add_fusion(fusion)
		
		self.assertEqual(len(mf_a), 2)
		
		#F[c] + F[d] = MF(C, D)
		mf_b = MergedFusion()
		for fusion in experiment_b:
			mf_b.add_fusion(fusion)
		
		self.assertEqual(len(mf_b), 2)
		
		self.assertEqual(str(mf_b), "--- MergedFusion container of size 2 ---\nFusion 't2_431' (from dataset 'test2'): chr22:15465000(-)<-chr22:41929200(+)\n\nFusion 't2_223' (from dataset 'test2'): chr11:524500(-)<-chr11:62910000(+)\n----------------------------------------\n")
		
		## Following functionality is deprecated for now:
		## MF(A , B) + MF(C , D) = MF(A, B, C, D)
		#mf_a.merge(mf_b)
		#del(mf_b)
		#self.assertEqual(len(mf_a), 4)




def main():
	unittest.main()

if __name__ == '__main__':
	main()
