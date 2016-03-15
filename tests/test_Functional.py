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

import unittest,subprocess,os,hashlib,sys,logging,subprocess
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",stream=sys.stdout)

class TestFusion(unittest.TestCase):
	def test_01(self):
		"""
		Functional test with test data
		"""
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command += ("bin/fuma \\\n"
						" --no-strand-specific-matching \\\n"
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
						" -m subset \\\n"
						" -f list \\\n"
						" -o test_Functional.test_01.output.txt "
					)
		
		os.system(command)
		
		# MD5 comparison:
		md5_input   = hashlib.md5(open('test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_01.output.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_Functional.test_01.output.txt')
	
	def test_02(self):
		#command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command = ["bin/fuma", \
					"--no-strand-specific-matching", \
					"-a","hg19:tests/data/refseq_hg19.bed", \
					"-s",
						"soapfuse-final-gene:soapfuse-final-gene:tests/data/test_Readers.TestReadSOAPFuseGenes.test_01.txt", \
						"soapfuse-final-transcript:soapfuse-final-transcript:tests/data/test_Readers.TestReadSOAPFuseTranscripts.test_01.txt", \
						"ericscript:ericscript:tests/data/test_Readers.TestReadEricScriptResultsTotal.test_01.txt", \
						"jaffa:jaffa:tests/data/test_Readers.TestReadJaffaResults.test_01.txt", \
					"-l", \
						"soapfuse-final-gene:hg19", \
						"soapfuse-final-transcript:hg19", \
						"ericscript:hg19", \
						"jaffa:hg19", \
					"-m","subset", \
					"-f","list", \
					"-o","-"]
		
		self.assertEqual(subprocess.call(command) , 0)# Ensure error code is 0 - no exceptions have been thrown
	
	def test_Edgren_hg19_summary(self):
		"""
		Functional test with test Edgren data (comparison to all genes on hg19)
		"""
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command += ("bin/fuma \\\n"
						" --no-strand-specific-matching \\\n"
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
						" -m subset \\\n"
						" -f summary \\\n"
						" -o test_Functional.test_Edgren_hg19.output.summary.txt "
					)
		
		os.system(command)
		
		# MD5 comparison:
		md5_input   = hashlib.md5(open('test_Functional.test_Edgren_hg19.output.summary.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_Edgren_hg19.output.summary.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_Functional.test_Edgren_hg19.output.summary.txt')
	
	def test_Edgren_hg19_list(self):
		"""
		Functional test with test Edgren data (comparison to all genes on hg19)
		"""
		
		command = "export PYTHONPATH=$PYTHONPATH\":fuma:../fuma\" ;\n\n"	# ensure the fuma lib is accessible for testing (also without installation)
		command += ("bin/fuma \\\n"
						" --no-strand-specific-matching \\\n"
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
						" -m subset \\\n"
						" -f list \\\n"
						" -o test_Functional.test_Edgren_hg19.output.list.txt"
					)
		
		os.system(command)
		
		# MD5 comparison:
		md5_input   = hashlib.md5(open('test_Functional.test_Edgren_hg19.output.list.txt', 'rb').read()).hexdigest()
		md5_confirm = hashlib.md5(open('tests/data/test_Functional.test_Edgren_hg19.output.list.txt', 'rb').read()).hexdigest()
		
		validation_1 = (md5_input != '')
		validation_2 = (md5_input == md5_confirm)
		
		self.assertNotEqual(md5_input , '')
		self.assertNotEqual(md5_confirm , '')
		self.assertEqual(md5_input , md5_confirm)
		
		if(validation_1 and validation_2):
			os.remove('test_Functional.test_Edgren_hg19.output.list.txt')


def main():
	unittest.main()

if __name__ == '__main__':
	main()
