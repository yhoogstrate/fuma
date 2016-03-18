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


import logging,fuma,math


class Triangle:
	logger = logging.getLogger("FuMa::ComparisonMatrix")
	
	def __init__(self,width=1,default=None):
		self.width = width
		self.data = [None] * len(self)
	
	def get(self,x,y):
		if not isinstance(arg_experiment, FusionDetectionExperiment):
			raise Exception("MergedFusion objects can only be expanded with Fusion objects")
		else:
			if arg_experiment in self.experiments:
				raise Exception("MergedFusion is updated with one that it already contains")
			else:
				self.experiments.append(arg_experiment)
	
	def overlay_fusions(self):
		triangle = self.create_triangle(len(self))
		
		#print list(self.experiments)[0].name
	
	def __len__(self):
		# n^2 - ((n-1) * n) / 2 = n * (n+1) / 2
		return (self.width * (self.width + 1)) / 2
	
def fsx(n):
	return n*(n+1)/2

def get_data_position(x,y):
	
	#m = 3
	#n = 9
	
	
	#n - math.factorial(m-y) + x
	
	
	loss_previous_rows = n - fsx(m-y)
	
	return loss_previous_rows
	
	"""
	for m = 3
	0 = 0
	1 = 3
	2 = 5
	#3 = 6
	
	for m = 4
	y=0 = 0       = 10 - (1 + 2 + 3 + 4) = fs(m-y)
	y=1 = 4       = 10 - (1 + 2 + 3)
	y=2 = 4+3     = 10 - (1 + 2)
	y=3 = 4+3+2   = 10 - (1)
	#4 = 4+3+2+1 = 10 - 
	"""
	bottom = ((r - 1) * r) / 2
	
	print r,bottom,row
	
	return n - bottom - row - 1
	
	def __str__(self):
		out = ""
		for i in range(self.width):
			out += "|"
			for j in range(i,self.width):
				out += "*|"
			out += "\n"
		
		print self.data
		
		return out
