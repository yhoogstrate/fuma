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

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.Readers import ReadDefuse
from fuma.Readers import ReadFusionMap
from fuma.ParseBED import ParseBED
from fuma.ComparisonTriangle import ComparisonTriangle
from fuma.CLI import CLI



def match_files_unsorted(filename_1, filename_2):
	fh1 = open(filename_1,"r")
	content1 = fh1.read().split("\n")
	fh1.close()
	
	fh2 = open(filename_2,"r")
	content2 = fh2.read().split("\n")
	fh2.close()
	
	content1.sort()
	content2.sort()
	
	if len(content1) == len(content2):
		for i in range(len(content1)):
			if content1[i] != content2[i]:
				return False
		
		return True
	else:
		return False


class TestComparisonTriangle(unittest.TestCase):
	def test_01(self):
		args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o','test_ComparisonTriangle.test_01.output.txt'])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		self.assertEqual(overlap.map_i_to_exp_id(0), 0)
		self.assertEqual(overlap.map_i_to_exp_id(1), 1)
		self.assertEqual(overlap.map_i_to_exp_id(2), 1)
		self.assertEqual(overlap.map_i_to_exp_id(3), 2)
		self.assertEqual(overlap.map_i_to_exp_id(4), 2)
		self.assertEqual(overlap.map_i_to_exp_id(5), 2)
		self.assertEqual(overlap.map_i_to_exp_id(6), 3)
		self.assertEqual(overlap.map_i_to_exp_id(7), 3)
		
		
		overlap.overlay_fusions()
		
				# MD5 comparison:
		md5_input   = hashlib.md5(open('test_ComparisonTriangle.test_01.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_ComparisonTriangle.test_01.output.txt')
		
		
	def test_02(self):
		args = CLI(['-m','egm','--no-strand-specific-matching','-s','','-o','test_ComparisonTriangle.test_02.output.txt'])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		self.assertEqual(overlap.map_i_to_exp_id(0), 0)
		self.assertEqual(overlap.map_i_to_exp_id(1), 1)
		self.assertEqual(overlap.map_i_to_exp_id(2), 1)
		self.assertEqual(overlap.map_i_to_exp_id(3), 2)
		self.assertEqual(overlap.map_i_to_exp_id(4), 2)
		self.assertEqual(overlap.map_i_to_exp_id(5), 2)
		self.assertEqual(overlap.map_i_to_exp_id(6), 3)
		self.assertEqual(overlap.map_i_to_exp_id(7), 3)
		
		
		overlap.overlay_fusions()
		
				# MD5 comparison:
		md5_input   = hashlib.md5(open('test_ComparisonTriangle.test_02.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_ComparisonTriangle.test_02.output.txt')
		
		
	def test_03(self):
		args = CLI(['-m','overlap','--no-strand-specific-matching','-s','','-o','test_ComparisonTriangle.test_03.output.txt'])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","test1")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","test2")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","test3")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","test4")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/refseq_hg19.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		self.assertEqual(len(overlap), 4)
		
		self.assertEqual(overlap.map_i_to_exp_id(0), 0)
		self.assertEqual(overlap.map_i_to_exp_id(1), 1)
		self.assertEqual(overlap.map_i_to_exp_id(2), 1)
		self.assertEqual(overlap.map_i_to_exp_id(3), 2)
		self.assertEqual(overlap.map_i_to_exp_id(4), 2)
		self.assertEqual(overlap.map_i_to_exp_id(5), 2)
		self.assertEqual(overlap.map_i_to_exp_id(6), 3)
		self.assertEqual(overlap.map_i_to_exp_id(7), 3)
		
		
		overlap.overlay_fusions()
		
		# MD5 comparison:
		md5_input   = hashlib.md5(open('test_ComparisonTriangle.test_03.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_ComparisonTriangle.test_03.output.txt')
	
	def test_04(self):
		"""
		Functional test with test Edgren data (comparison to all genes on hg19)
		"""
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o','test_ComparisonTriangle.test_04.output.txt'])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_Edgren_hg19.ChimeraScan.txt","chimerascan")
		experiment_b = ReadDefuse("tests/data/test_Functional.test_Edgren_hg19.Defuse.txt","defuse")
		experiment_c = ReadFusionMap("tests/data/test_Functional.test_Edgren_hg19.FusionMap.txt","fusion-map")
		experiment_d = ReadFusionMap("tests/data/test_Functional.test_Edgren_hg19.TruePositives.txt","edgren_tp")
		
		genes = ParseBED("tests/data/refseq_genes_hg19.bed","hg19",200000)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_a)
		overlap.add_experiment(experiment_b)
		overlap.add_experiment(experiment_c)
		overlap.add_experiment(experiment_d)
		
		
		overlap.overlay_fusions()
		
		## Order may have changed with respect to earlier releases, but the same fusion genes are reported
		files_identical = match_files_unsorted('test_ComparisonTriangle.test_04.output.txt','tests/data/test_Functional.test_Edgren_hg19.output.list.txt')
		self.assertTrue(files_identical)
		
		if(not files_identical):
			os.remove('test_ComparisonTriangle.test_04.output.txt')



def main():
	unittest.main()

if __name__ == '__main__':
	main()
