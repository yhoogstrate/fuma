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

from Readers import *

from ParseBED import ParseBED
from FusionDetectionExperiment import FusionDetectionExperiment
from Triangle import Triangle


import os.path,sys,itertools


class ComparisonMatrix:
	logger = logging.getLogger("FuMa::ComparisonMatrix")
	
	def __init__(self,args):
		self.experiments = []
		self.args = args
	
	def add_experiment(self,arg_experiment):
		if not isinstance(arg_experiment, FusionDetectionExperiment):
			raise Exception("MergedFusion objects can only be expanded with Fusion objects")
		else:
			if arg_experiment in self.experiments:
				raise Exception("MergedFusion is updated with one that it already contains")
			else:
				self.experiments.append(arg_experiment)
	
	def overlay_fusions(self):
		triangle = Triangle(len(self))
		
		for experiment in self.experiments:
			for fusion in 
	
	def __len__(self):
		return len(self.experiments)
