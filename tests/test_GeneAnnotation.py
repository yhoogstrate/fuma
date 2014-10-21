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

from fuma.Gene import Gene
from fuma.GeneAnnotation import GeneAnnotation

class TestGeneAnnotation(unittest.TestCase):
	def test_01(self):
		genes = GeneAnnotation("hg18")
		
		gene_01 = Gene("ucsc.1")
		gene_02 = Gene("ucsc.2")
		gene_03 = Gene("ucsc.3")
		gene_04 = Gene("ucsc.4")
		
		genes.add_annotation(gene_01,"chr3",10,15)
		genes.add_annotation(gene_02,"chr3",11,16)
		genes.add_annotation(gene_03,"chr3",12,18)
		genes.add_annotation(gene_04,"chr3",12,20)
		
		self.assertEqual(len(genes), 4)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
