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

import logging,os



from Fusion import Fusion
from FusionDetectionExperiment import FusionDetectionExperiment




def euclidean(x,y):
	sumSq = 0.0
	
	#add up the squared differences
	for i in range(len(x)):
		sumSq += (x[i]-y[i])**2
		
	#take the square root of the result
	return (sumSq**0.5)



class CompareFusionsDistance:
	def __init__(self,dataset_1,dataset_2):
		self.dataset_1 = dataset_1
		self.dataset_2 = dataset_2
		self.matches = None
	
	def find_all_distances(self):
		"""
		n1 * n2 comparisons -> n1 * n2 results
		"""
		print " - Comparing: "+self.dataset_1.name+" - "+self.dataset_2.name
		
		self.matches = []
		
		for chromosome_1 in self.dataset_1.get_fusions():
			for fusion_1 in chromosome_1["fusions"]:
				
				for chromosome_2 in self.dataset_2.get_fusions():
					if(chromosome_1["name"] == chromosome_2["name"]):
						for fusion_2 in chromosome_2["fusions"]:
							self.matches.append(self.find_distance(fusion_1,fusion_2))
		return self.matches
	
	def find_minimum_distances(self):
		"""
		n1 * n2 comparisons -> n1
		"""
		print " - Comparing: "+self.dataset_1.name+" - "+self.dataset_2.name
		
		self.matches = []
		
		for chromosome_1 in self.dataset_1.get_fusions():
			for fusion_1 in chromosome_1["fusions"]:
				
				distance = {"fusion_1":fusion_1,"fusion_2":None,"city block distance":"inf","eucledian distance":"inf"}
				
				for chromosome_2 in self.dataset_2.get_fusions():
					if(chromosome_1["name"] == chromosome_2["name"]):
						for fusion_2 in chromosome_2["fusions"]:
							tmp_distance = self.find_distance(fusion_1,fusion_2)
							if(tmp_distance["eucledian distance"] != "inf"):
								if(distance["eucledian distance"] == "inf" or (tmp_distance["eucledian distance"] < distance["eucledian distance"])):
									distance = tmp_distance
				
				self.matches.append(distance)
		return self.matches
	
	def find_distance(self,fusion_1,fusion_2):
		if((fusion_1.get_left_chromosome(False) != fusion_2.get_left_chromosome(False)) or (fusion_1.get_right_chromosome(False) != fusion_2.get_right_chromosome(False))):
			return {"fusion_1":fusion_1,"fusion_2":fusion_2,"city block distance":"inf","eucledian distance":"inf"}
		else:
			left_distance = abs(fusion_1.get_left_break_position() - fusion_2.get_left_break_position())
			right_distance = abs(fusion_1.get_right_break_position() - fusion_2.get_right_break_position())
			
			city_block_distance = left_distance + right_distance
			eucledian_distance = euclidean([fusion_1.get_left_break_position(),fusion_1.get_right_break_position()],[fusion_2.get_left_break_position(),fusion_2.get_right_break_position()])
			
			return {"fusion_1":fusion_1,"fusion_2":fusion_2,"city block distance":city_block_distance,"eucledian distance":eucledian_distance}
	
	def export(self,filename_prefix="",join="_vs._",suffix=".txt"):
		if(self.matches):
			filename = filename_prefix+self.dataset_1.name+join+os.path.basename(self.dataset_2.name)+suffix
			print "exporting: "+os.path.basename(filename)
			fh = open(filename,"w")
			fh.write("left_fusion["+self.dataset_1.name+"]\t")
			fh.write("right_fusion["+self.dataset_2.name+"]\t")
			fh.write("eucledian distance\t")
			fh.write("city block distance\n")
			
			keys = {}
			for distance in self.matches:
				#print distance
				
				
				if(not keys.has_key(distance["eucledian distance"])):
					keys[distance["eucledian distance"]] = []
				keys[distance["eucledian distance"]].append(distance)
			
			for key in sorted(keys.keys()):
				for distance in keys[key]:
					if distance["fusion_2"]:
						fh.write(str(distance["fusion_1"].get_left_chromosome())+":")
						fh.write(str(distance["fusion_1"].get_left_break_position())+"-")
						fh.write(str(distance["fusion_1"].get_right_chromosome())+":")
						fh.write(str(distance["fusion_1"].get_right_break_position())+"\t")
						
						fh.write(str(distance["fusion_2"].get_left_chromosome())+":")
						fh.write(str(distance["fusion_2"].get_left_break_position())+"-")
						fh.write(str(distance["fusion_2"].get_right_chromosome())+":")
						fh.write(str(distance["fusion_2"].get_right_break_position())+"\t")
						
						fh.write(str(distance["eucledian distance"])+"\t")
						fh.write(str(distance["city block distance"])+"\n")
					else:
						fh.write(str(distance["fusion_1"].get_left_chromosome())+":")
						fh.write(str(distance["fusion_1"].get_left_break_position())+"-")
						fh.write(str(distance["fusion_1"].get_right_chromosome())+":")
						fh.write(str(distance["fusion_1"].get_right_break_position())+"\t")
						
						fh.write("\t\t\n")
			fh.close()
