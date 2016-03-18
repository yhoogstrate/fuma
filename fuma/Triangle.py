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
		return self.data[self.get_data_position(x,y)]
	
	def set(self,x,y,value):
		self.data[self.get_data_position(x,y)] = value
	
	def __len__(self):
		# n^2 - ((n-1) * n) / 2 = n * (n+1) / 2
		return (self.width * (self.width + 1)) / 2
	
	def get_data_position(self,x,y):
		def fsx(n):
			return n*(n+1)/2
		return len(self.data) - fsx(self.width-y) + x
	
	def __str__(self):
		out = ""
		for i in range(self.width):
			out += "|"
			for j in range(i,self.width):
				out += "*|"
			out += "\n"
		
		print self.data
		
		return out
