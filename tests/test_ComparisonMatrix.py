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
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.Readers import ReadChimeraScanAbsoluteBEDPE
from fuma.ParseBED import ParseBED
from fuma.ComparisonMatrix import ComparisonMatrix
from fuma.CLI import CLI

class TestComparisonMatrix(unittest.TestCase):
	def test_01(self):
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","TestExperimentB")
		experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","TestExperimentC")
		experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","TestExperimentD")
		
		self.assertEqual(len(experiment_a), 2)
		self.assertEqual(len(experiment_b), 2)
		self.assertEqual(len(experiment_c), 3)
		self.assertEqual(len(experiment_d), 3)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg19",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		experiment_c.annotate_genes(genes)
		experiment_d.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		experiment_c.remove_duplicates(args)
		experiment_d.remove_duplicates(args)
		
		overlap = ComparisonMatrix(args)
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
		
		
	#def test_02(self):
		#args = CLI(['-m','egm','--no-strand-specific-matching','-s',''])
		
		#experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","TestExperimentA")
		#experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","TestExperimentB")
		#experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","TestExperimentC")
		#experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","TestExperimentD")
		
		#self.assertEqual(len(experiment_a), 2)
		#self.assertEqual(len(experiment_b), 2)
		#self.assertEqual(len(experiment_c), 3)
		#self.assertEqual(len(experiment_d), 3)
		
		#genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg19",200000)
		
		#self.assertEqual(len(genes), 47790)
		
		#experiment_a.annotate_genes(genes)
		#experiment_b.annotate_genes(genes)
		#experiment_c.annotate_genes(genes)
		#experiment_d.annotate_genes(genes)
		
		#experiment_a.remove_duplicates(args)
		#experiment_b.remove_duplicates(args)
		#experiment_c.remove_duplicates(args)
		#experiment_d.remove_duplicates(args)
		
		#overlap = ComparisonMatrix(args)
		#overlap.add_experiment(experiment_a)
		#overlap.add_experiment(experiment_b)
		#overlap.add_experiment(experiment_c)
		#overlap.add_experiment(experiment_d)
		
		#self.assertLessEqual(len(overlap), 2+2+3+3)
		
		#overlap.overlay_fusions()
		
		
	#def test_03(self):
		#args = CLI(['-m','overlap','--no-strand-specific-matching','-s',''])
		
		#experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_01.bedpe","TestExperimentA")
		#experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_02.bedpe","TestExperimentB")
		#experiment_c = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_03.bedpe","TestExperimentC")
		#experiment_d = ReadChimeraScanAbsoluteBEDPE("tests/data/test_Functional.test_01.Example_04.bedpe","TestExperimentD")
		
		#self.assertEqual(len(experiment_a), 2)
		#self.assertEqual(len(experiment_b), 2)
		#self.assertEqual(len(experiment_c), 3)
		#self.assertEqual(len(experiment_d), 3)
		
		#genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg19",200000)
		
		#self.assertEqual(len(genes), 47790)
		
		#experiment_a.annotate_genes(genes)
		#experiment_b.annotate_genes(genes)
		#experiment_c.annotate_genes(genes)
		#experiment_d.annotate_genes(genes)
		
		#experiment_a.remove_duplicates(args)
		#experiment_b.remove_duplicates(args)
		#experiment_c.remove_duplicates(args)
		#experiment_d.remove_duplicates(args)
		
		#overlap = ComparisonMatrix(args)
		#overlap.add_experiment(experiment_a)
		#overlap.add_experiment(experiment_b)
		#overlap.add_experiment(experiment_c)
		#overlap.add_experiment(experiment_d)
		
		#self.assertLessEqual(len(overlap), 2+2+3+3)
		
		#overlap.overlay_fusions()


def main():
	unittest.main()

if __name__ == '__main__':
	main()
