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

import unittest,subprocess,os,hashlib

class TestFusion(unittest.TestCase):
	def test_01(self):
		"""
		Functional test with test data
		"""
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing also without installation
		command += ("bin/fuma \\\n"
						" -a hg19:tests/data/refseq_hg19.bed \\\n"
						" -s \\\n"
						"   test1:chimerascan:tests/data/test_Functional.test_01.Example_01.bedpe \\\n"
						"   test2:chimerascan:tests/data/test_Functional.test_01.Example_02.bedpe \\\n"
						"   test3:chimerascan:tests/data/test_Functional.test_01.Example_03.bedpe \\\n"
						"   test4:chimerascan:tests/data/test_Functional.test_01.Example_04.bedpe \\\n"
						" -l \\\n"
						"    test1:hg19 \\\n"
						"    test2:hg19 \\\n"
						"    test3:hg19 \\\n"
						"    test4:hg19 \\\n"
						" -f list \\\n"
						" -o test_Functional.test_01.output.txt "
					)
		
		os.system(command)
		
		md5_input   = hashlib.md5(open('test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		
		self.assertNotEqual(md5_input, '')
		self.assertNotEqual(md5_confirm, '')
		
		self.assertEqual(md5_input , md5_confirm)
		
		if(md5_input == md5_confirm):
			os.remove('test_Functional.test_01.output.txt')
	
	
	def test_Edgren_hg19(self):
		"""
		Functional test with test Edgren data (comparison to all genes on hg19)
		"""
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing also without installation
		command += ("bin/fuma \\\n"
						" -a hg19:tests/data/refseq_genes_hg19.bed \\\n"
						" -s \\\n"
						"   chimerascan:chimerascan:tests/data/test_Functional.test_Edgren_hg19.ChimeraScan.txt \\\n"
						"             defuse:defuse:tests/data/test_Functional.test_Edgren_hg19.Defuse.txt \\\n"
						"      fusion-map:fusionmap:tests/data/test_Functional.test_Edgren_hg19.FusionMap.txt \\\n"
						"       edgren_tp:fusionmap:tests/data/test_Functional.test_Edgren_hg19.TruePositives.txt \\\n"
						" -l \\\n"
						"    chimerascan:hg19 \\\n"
						"         defuse:hg19 \\\n"
						"     fusion-map:hg19 \\\n"
						"      edgren_tp:hg19 \\\n"
						" -f summary \\\n"
						" -o test_Functional.test_Edgren_hg19.output.txt "
					)
		
		os.system(command)
		
		md5_input   = hashlib.md5(open('test_Functional.test_Edgren_hg19.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_Edgren_hg19.output.txt', 'rb').read()).hexdigest()
		
		self.assertNotEqual(md5_input, '')
		self.assertNotEqual(md5_confirm, '')
		
		self.assertEqual(md5_input , md5_confirm)
		
		if(md5_input == md5_confirm):
			os.remove('test_Functional.test_Edgren_hg19.output.txt')


def main():
	unittest.main()

if __name__ == '__main__':
	main()



