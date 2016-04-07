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


from fuma.CLI import CLI
from fuma.Readers import ReadRNASTARFusionFinal
from fuma.ParseBED import ParseBED
from fuma.ComparisonTriangle import ComparisonTriangle


class TestListOutput(unittest.TestCase):
	def test_01(self):
		output_file = 'test_ListOutput.test_01.output.txt'
		validation_file = 'tests/data/test_ListOutput.test_01.output.txt'
		
		args = CLI(['-m','subset','--no-strand-specific-matching','-s','','-o',output_file])
		
		experiment_41 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-041.txt","test-41")
		experiment_42 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-042.txt","test-42")
		experiment_43 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-043.txt","test-43")
		experiment_44 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-044.txt","test-44")
		experiment_45 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-045.txt","test-45")
		experiment_47 = ReadRNASTARFusionFinal("tests/high_dim/7046-004-047.txt","test-47")
		
		genes = ParseBED("tests/high_dim/hg19_subset.bed","hg19",200000)
		
		experiment_41.annotate_genes(genes)
		experiment_42.annotate_genes(genes)
		experiment_43.annotate_genes(genes)
		experiment_44.annotate_genes(genes)
		experiment_45.annotate_genes(genes)
		experiment_47.annotate_genes(genes)
		
		experiment_41.remove_duplicates(args)
		experiment_42.remove_duplicates(args)
		experiment_43.remove_duplicates(args)
		experiment_44.remove_duplicates(args)
		experiment_45.remove_duplicates(args)
		experiment_47.remove_duplicates(args)
		
		overlap = ComparisonTriangle(args)
		overlap.add_experiment(experiment_41)
		overlap.add_experiment(experiment_42)
		overlap.add_experiment(experiment_43)
		overlap.add_experiment(experiment_44)
		overlap.add_experiment(experiment_45)
		overlap.add_experiment(experiment_47)
		overlap.overlay_fusions()
		
		# mix up order
		#self.assertEqual( 1,1 )

def main():
	unittest.main()

if __name__ == '__main__':
	main()
