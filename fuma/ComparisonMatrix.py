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
		#triangle = Triangle(len(self))
		fusions = []
		
		for experiment in self.experiments:
			for fusion in experiment:
				fusions.append(fusion)
		
		n = len(fusions)
		
		"""
		iter1: 
		
		0,7 | 1,7 2,7 | 3,7 4,7 5,7 | 6,7 7,7
		
		iter2:
		0,6 | 1,6 2,6 | 3,6 4,6 5,6 | 6,6
		
		"""
		for y in range(len(fusions)-1,0,-1):
			fusions_y = fusions[y]
			for x in range(y+1):
				if self.map_i_to_exp_id(x) != self.map_i_to_exp_id(y):# If they do not belong to the same dataset - i.e. no duplication removal
					fusion_x = fusions[x]
					
					## do actual comarison
			
		
		"""
		should be in unit test:
		
		print "0="+str(self.map_i_to_exp_id(0))+" =? 0"
		print
		print "1="+str(self.map_i_to_exp_id(1))+" =? 1"
		print
		print "2="+str(self.map_i_to_exp_id(2))+" =? 1"
		print
		print "3="+str(self.map_i_to_exp_id(3))+" =? 2"
		print
		print "4="+str(self.map_i_to_exp_id(4))+" =? 2"
		print
		print "5="+str(self.map_i_to_exp_id(5))+" =? 2"
		print
		print "6="+str(self.map_i_to_exp_id(6))+" =? 3"
		print
		print "7="+str(self.map_i_to_exp_id(7))+" =? 3"
		print
		"""
	
	def map_i_to_exp_id(self,i):
		# 1 2 3 2
		cumulative = 0
		j = -1
		while i >= cumulative:
			j += 1
			cumulative += len(self.experiments[j])
		return j
	
	def __len__(self):
		return len(self.experiments)
