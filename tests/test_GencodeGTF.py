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

import unittest,logging,sys,os
logging.basicConfig(level=logging.INFO,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

from fuma.ParseBED import ParseBED

class TestParseBED(unittest.TestCase):
	def test_01(self):
		inputfile = "tests/data/gencode_hg19.subset.gtf"
		outputfile = "tests/data/gencode_hg19.subset.bed"
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command += ("bin/fuma-gencode-gtf-to-bed\\\n"
					"   "+inputfile
					)
		
		result = os.popen(command).read()
		validation = open(outputfile,"r").read()
		
		self.assertEqual(result, validation)

	def test_02(self):
		inputfile = "tests/data/gencode_hg19.subset.gtf"
		outputfile = "tests/data/gencode_hg19.subset.sorted.bed"
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command += ("bin/fuma-gencode-gtf-to-bed \\\n"
					"   "+inputfile+" | sort -k1,1V -k2,2g -k3,3g "
					)

		
		result = os.popen(command).read()
		validation = open(outputfile,"r").read()
		
		self.assertEqual(result, validation)

def main():
	unittest.main()

if __name__ == '__main__':
	main()
