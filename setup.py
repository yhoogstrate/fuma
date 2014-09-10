	#!/usr/bin/env python

	"""[License: GPL3]
	 
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

	import FuMa

	from distutils.core import setup

	setup(name='fuma',
		  version=fuma.__version__,
		  description='Fusion Matcher',
		  author=fuma.__author__,
		  url='https://github.com/yhoogstrate/fuma',
		  scripts=["bin/fuma"],
		  packages=['fuma'],
		 )
