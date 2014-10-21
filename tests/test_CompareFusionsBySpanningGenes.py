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

class TestCompareFusionsBySpanningGenes(unittest.TestCase):
	def test_01(self):
		experiment_a = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentA")
		experiment_b = ReadChimeraScanAbsoluteBEDPE("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bedpe","TestExperimentB")
		
		self.assertEqual(len(experiment_a), 690)
		self.assertEqual(len(experiment_b), 690)
		
		genes = ParseBED("tests/data/test_CompareFusionsBySpanningGenes.TestCompareFusionsBySpanningGenes.test_01.bed","hg18")
		
		self.assertEqual(len(genes), 47790)
		
		experiment_a.annotate_genes(genes)
		experiment_b.annotate_genes(genes)
		
		experiment_a.remove_duplicates("by-gene-names")
		experiment_b.remove_duplicates("by-gene-names")
		
		overlap = CompareFusionsBySpanningGenes(experiment_a,experiment_b)
		overlapping_fusions = overlap.find_overlap()
		
		self.assertLessEqual(len(overlapping_fusions), len(experiment_a))
		
		#@todo test what the difference is between the outcome of this version and the classical version(!) on the same datasets
		
		#self.assertEqual(i, 690)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
