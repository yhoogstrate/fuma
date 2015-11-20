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
from fuma.CompareFusionsBySpanningGenes import CompareFusionsBySpanningGenes
from fuma.CLI import CLI

class TestCompareFusionsBySpanningGenes(unittest.TestCase):
	def test_01(self):
		args = CLI(['-m','subset','--no-strand-specific-matching','-s',''])
		
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentB")
		
		self.assertEqual(len(experiment_a), 690)
		self.assertEqual(len(experiment_b), 690)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg18",200000)
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		
		experiment_a.remove_duplicates(args)
		experiment_b.remove_duplicates(args)
		
		overlap = CompareFusionsBySpanningGenes(experiment_a,experiment_b,args)
		overlapping_fusions = overlap.find_overlap()
		
		self.assertLessEqual(len(overlapping_fusions[0]), len(experiment_a))
	
	def test_02(self):
		args_a = CLI(['-m','subset','-s',''])
		args_b = CLI(['-m','subset','--strand-specific-matching','-s',''])
		
		## First test the matches if strand-specific-matching is disabled (all 4 fusions should be identical)
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02_a.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02_b.bedpe","TestExperimentB")
		
		self.assertEqual(len(experiment_a), 4)
		self.assertEqual(len(experiment_b), 4)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_02.bed","hg18",200000)
		
		self.assertEqual(len(genes), 8)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		
		## @todo -> remove duplicates should be done separately
		experiment_a.remove_duplicates(args_a)
		experiment_b.remove_duplicates(args_a)
		
		overlap = CompareFusionsBySpanningGenes(experiment_a,experiment_b,args_a)# No EGM, no strand-specific-matching
		overlapping_fusions = overlap.find_overlap()
		
		self.assertLessEqual(len(overlapping_fusions[0]), 4)
		
		## Second, test the matches if strand-specific-matching is disabled (only the first fusion should be identical)
		overlap = CompareFusionsBySpanningGenes(experiment_a,experiment_b,args_b)# No EGM, but strand-specific-matching
		overlapping_fusions = overlap.find_overlap()
		
		self.assertLessEqual(len(overlapping_fusions[0]), 1)


def main():
	unittest.main()

if __name__ == '__main__':
	main()
