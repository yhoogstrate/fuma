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

from fuma.ParseBED import ParseBED

class TestParseBED(unittest.TestCase):
	def test_01(self):
		filename = "tests/data/test_ParseBED.TestParseBED.test_01.bed"
		gene_annotation = ParseBED(filename,"hg1234")
		#gene_annotation.show_me()
		#self.failUnless(1 == 1)
		#self.failIf(1 != 1)
		#

def main():
	unittest.main()

if __name__ == '__main__':
	main()
