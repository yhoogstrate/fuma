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

import fuma

from setuptools import setup

setup(name='fuma',
		version=fuma.__version__,
		description='Fusion Matcher',
		long_description="FuMa (Fusion Matcher) matches predicted fusion events (both genomic and transcriptomic) according to chromosomal location and corresponding annotated genes. It is the organisation of the transcriptome (provided by the user) that forms the basis for FuMa to consider fusion genes to be identical or not. The provided gene annotation can be adjusted to define the biological question. For example, if it is desired to only consider fusion events that occur within exons, FuMa can be provided a list of such regions instead of entire genes.",
		author=fuma.__author__,
		maintainer=fuma.__author__,
		license=fuma.__license__,
		url=fuma.__homepage__,
		scripts=["bin/fuma","bin/defuse-clusters-to-CG",'bin/chimerascan-exclude-transcriptome-events',"bin/fusioncatcher-to-CG","bin/chimerascan-relative-bedpe-to-CG","bin/fuma-list-to-boolean-list","bin/fuma-gencode-gtf-to-bed"],
		packages=['fuma'],
		test_suite="tests",
		platforms=['any'],
		setup_requires=['numpy'],
		install_requires=['numpy','HTSeq >= 0.6.1','nose'],
		classifiers=[
			'Environment :: Console',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
			'Operating System :: OS Independent',
			'Topic :: Scientific/Engineering',
			'Topic :: Scientific/Engineering :: Bio-Informatics'
			],
	)
