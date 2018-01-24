#!/usr/bin/env python 
"""
Copyright (C) 2009 Nick Drobchenko, nick@cnc-club.ru
based on gcode.py (C) 2007 hugomatic... 
based on addnodes.py (C) 2005,2007 Aaron Spike, aaron@ekips.org
based on dots.py (C) 2005 Aaron Spike, aaron@ekips.org
based on interp.py (C) 2005 Aaron Spike, aaron@ekips.org
based on bezmisc.py (C) 2005 Aaron Spike, aaron@ekips.org
based on cubicsuperpath.py (C) 2005 Aaron Spike, aaron@ekips.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

###
###		Unordinary gears v. 1.1
###


import inkex, simplestyle, simplepath
import cubicsuperpath, simpletransform, bezmisc

import os
import math
import bezmisc
import re
import copy
import sys
import time
import cmath
import numpy
import codecs
import random
import gettext
import time
_ = gettext.gettext


optimize_distance_max_iters = 25
optimize_distance_tolerance = 0.00001

 ### Check if inkex has errormsg (0.46 version doesnot have one.) Could be removed later.
if "errormsg" not in dir(inkex):
	inkex.errormsg = lambda msg: sys.stderr.write((unicode(msg) + "\n").encode("UTF-8"))

def csp_at_t(sp1,sp2,t):
	ax,bx,cx,dx = sp1[1][0], sp1[2][0], sp2[0][0], sp2[1][0]
	ay,by,cy,dy = sp1[1][1], sp1[2][1], sp2[0][1], sp2[1][1]

	x1, y1 = ax+(bx-ax)*t, ay+(by-ay)*t	
	x2, y2 = bx+(cx-bx)*t, by+(cy-by)*t	
	x3, y3 = cx+(dx-cx)*t, cy+(dy-cy)*t	
	
	x4,y4 = x1+(x2-x1)*t, y1+(y2-y1)*t 
	x5,y5 = x2+(x3-x2)*t, y2+(y3-y2)*t 
	
	x,y = x4+(x5-x4)*t, y4+(y5-y4)*t 
	return [x,y]


def csp_true_bounds(csp) :
	# Finds minx,miny,maxx,maxy of the csp and return their (x,y,i,j,t) 
	minx = [float("inf"),  0, 0, 0]																																																																																																								
	maxx = [float("-inf"), 0, 0, 0]
	miny = [float("inf"),  0, 0, 0]
	maxy = [float("-inf"), 0, 0, 0]
	for i in range(len(csp)):
		for j in range(1,len(csp[i])):
			ax,ay,bx,by,cx,cy,x0,y0 = bezmisc.bezierparameterize((csp[i][j-1][1],csp[i][j-1][2],csp[i][j][0],csp[i][j][1]))
			roots = cubic_solver(0, 3*ax, 2*bx, cx)	 + [0,1]
			for root in roots :
				if type(root) is complex and abs(root.imag)<1e-10:
					root = root.real
				if type(root) is not complex and 0<=root<=1:
					y = ay*(root**3)+by*(root**2)+cy*root+y0  
					x = ax*(root**3)+bx*(root**2)+cx*root+x0  
					maxx = max([x,y,i,j,root],maxx)
					minx = min([x,y,i,j,root],minx)

			roots = cubic_solver(0, 3*ay, 2*by, cy)	 + [0,1]
			for root in roots :
				if type(root) is complex and root.imag==0:
					root = root.real
				if type(root) is not complex and 0<=root<=1:
					y = ay*(root**3)+by*(root**2)+cy*root+y0  
					x = ax*(root**3)+bx*(root**2)+cx*root+x0  
					maxy = max([y,x,i,j,root],maxy)
					miny = min([y,x,i,j,root],miny)
	maxy[0],maxy[1] = maxy[1],maxy[0]
	miny[0],miny[1] = miny[1],miny[0]

	return minx,miny,maxx,maxy

def csp_segment_to_bez(sp1,sp2) :
	return sp1[1:]+sp2[:2]

def csp_parameterize(sp1,sp2):
	return bezmisc.bezierparameterize(csp_segment_to_bez(sp1,sp2))

def csp_to_point_distance(csp, p, dist_bounds = [0,1e100], tolerance=.001) :
	min_dist = [1e100,0,0,0]
	for j in range(len(csp)) :
		for i in range(1,len(csp[j])) :
			d = csp_seg_to_point_distance(csp[j][i-1],csp[j][i],p,sample_points = 5, tolerance = .001)
			if d[0] < dist_bounds[0] : 
#				draw_pointer( list(csp_at_t(subpath[dist[2]-1],subpath[dist[2]],dist[3]))
#					+list(csp_at_t(csp[dist[4]][dist[5]-1],csp[dist[4]][dist[5]],dist[6])),"red","line", comment = math.sqrt(dist[0]))
				return [d[0],j,i,d[1]]
			else : 
				if d[0] < min_dist[0] : min_dist = [d[0],j,i,d[1]]
	return min_dist
			
def csp_seg_to_point_distance(sp1,sp2,p,sample_points = 5, tolerance = .001) :
	ax,ay,bx,by,cx,cy,dx,dy = csp_parameterize(sp1,sp2)
	dx, dy = dx-p[0], dy-p[1]
	if sample_points < 2 : sample_points = 2
	d = min( [(p[0]-sp1[1][0])**2 + (p[1]-sp1[1][1])**2,0.], [(p[0]-sp2[1][0])**2 + (p[1]-sp2[1][1])**2,1.]	)	
	for k in range(sample_points) :
		t = float(k)/(sample_points-1)
		i = 0
		while i==0 or abs(f)>0.000001 and i<20 : 
			t2,t3 = t**2,t**3
			f = (ax*t3+bx*t2+cx*t+dx)*(3*ax*t2+2*bx*t+cx) + (ay*t3+by*t2+cy*t+dy)*(3*ay*t2+2*by*t+cy)
			df = (6*ax*t+2*bx)*(ax*t3+bx*t2+cx*t+dx) + (3*ax*t2+2*bx*t+cx)**2 + (6*ay*t+2*by)*(ay*t3+by*t2+cy*t+dy) + (3*ay*t2+2*by*t+cy)**2
			if df!=0 :
				t = t - f/df
			else :	
				break
			i += 1	
		if 0<=t<=1 : 
			p1 = csp_at_t(sp1,sp2,t)
			d1 = (p1[0]-p[0])**2 + (p1[1]-p[1])**2
			if d1 < d[0] :
				d = [d1,t]
	return d	



	
def cubic_solver(a,b,c,d):	
	if a!=0:
		#	Monics formula see http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
		a,b,c = (b/a, c/a, d/a)
		m = 2*a**3 - 9*a*b + 27*c
		k = a**2 - 3*b
		n = m**2 - 4*k**3
		w1 = -.5 + .5*cmath.sqrt(3)*1j
		w2 = -.5 - .5*cmath.sqrt(3)*1j
		if n>=0 :
			t = m+math.sqrt(n)
			m1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
			t = m-math.sqrt(n)
			n1 = pow(t/2,1./3) if t>=0 else -pow(-t/2,1./3)
		else :
			m1 = pow(complex((m+cmath.sqrt(n))/2),1./3)
			n1 = pow(complex((m-cmath.sqrt(n))/2),1./3)
		x1 = -1./3 * (a + m1 + n1)
		x2 = -1./3 * (a + w1*m1 + w2*n1)
		x3 = -1./3 * (a + w2*m1 + w1*n1)
		return [x1,x2,x3]
	elif b!=0:
		det = c**2-4*b*d
		if det>0 :
			return [(-c+math.sqrt(det))/(2*b),(-c-math.sqrt(det))/(2*b)]
		elif d == 0 :
			return [-c/(b*b)] 	
		else :
			return [(-c+cmath.sqrt(det))/(2*b),(-c-cmath.sqrt(det))/(2*b)]
	elif c!=0 :
		return [-d/c]
	else : return []	
	

class Unordinary_gears(inkex.Effect):
	
	def get_transforms(self,g):
		root = self.document.getroot()
		trans = []
		while (g!=root):
			if 'transform' in g.keys():
				t = g.get('transform')
				t = simpletransform.parseTransform(t)
				trans = simpletransform.composeTransform(t,trans) if trans != [] else t
			g=g.getparent()
		return trans 
	

	def apply_transforms(self,g,csp):
		trans = self.get_transforms(g)
		if trans != []:
			simpletransform.applyTransformToPath(trans, csp)
		return csp

	def __init__(self):
		inkex.Effect.__init__(self)
		self.OptionParser.add_option("", "--selected_puley_teeth",		action="store", type="int", 		dest="selected_puley_teeth", default="10",		help="Selected pulley teeth")
		self.OptionParser.add_option("", "--generated_puley_teeth",		action="store", type="int", 		dest="generated_puley_teeth", default="20",		help="Generated pulley teeth")
		self.OptionParser.add_option("", "--number_of_copies",			action="store", type="int", 		dest="number_of_copies", default="50",			help="Number of copies")
		self.OptionParser.add_option("", "--distance",					action="store", type="float", 		dest="distance", default="100",					help="Center to center distance")
		self.OptionParser.add_option("", "--units",						action="store", type="float", 		dest="units", default="90",						help="Units")
		self.OptionParser.add_option("", "--variable_speed",			action="store", type="inkbool", 	dest="variable_speed", default=False,			help="Turning speed depends on local radius")
		self.OptionParser.add_option("", "--optimize_distance",			action="store", type="inkbool", 	dest="optimize_distance", default=False,		help="Optimize center to center distance")

		self.OptionParser.add_option("", "--linear_distance",			action="store", type="float", 		dest="linear_distance", default=100.,			help="Center to shaft distance")
		self.OptionParser.add_option("", "--rev",						action="store", type="float", 		dest="rev", default=False,						help="Revolutions")

		self.OptionParser.add_option("", "--create_bearing_gear",		action="store", type="inkbool", 	dest="create_bearing_gear", default=False,		help="Create bearing gear")
		self.OptionParser.add_option("", "--number_of_bearings",		action="store", type="int",			dest="number_of_bearings", default=5,			help="Number of bearings")
		self.OptionParser.add_option("", "--bearings_d",				action="store", type="float", 		dest="bearings_d", default=15.,				help="Bearings diameter")
		self.OptionParser.add_option("", "--gear_r",					action="store", type="float", 		dest="gear_r", default=100.,					help="Gear's radius")
		self.OptionParser.add_option("", "--standatd_distance",			action="store", type="inkbool", 	dest="standatd_distance", default=False,		help="Use Gear R + Bearing R as distance")

		self.OptionParser.add_option("", "--active_tab",				action="store", type="string", 		dest="active_tab", default="",					help="")


	def error(self, s, type_= "Warning"):
		s = str(s)
		warnings = """
						Warning						
						"""
		errors = """
						Error 	
					"""
		if type_.lower() in re.split("[\s\n,\.]+", errors.lower()) :
			inkex.errormsg(s+"\n")		
			sys.exit()
		elif type_.lower() in re.split("[\s\n,\.]+", warnings.lower()) :
			inkex.errormsg(s+"\n")		
		else :
			inkex.errormsg(s)		
			sys.exit()

	def draw_pointer(self, x,color = "#f00", figure = "cross", comment = "", width = .1, group = None) :
		if group == None :
			group = self.document.getroot()
		if figure ==  "line" :
			s = ""
			for i in range(1,len(x)/2) :
				s+= " %s, %s " %(x[i*2],x[i*2+1])
			inkex.etree.SubElement( group, inkex.addNS('path','svg'), {"d": "M %s,%s L %s"%(x[0],x[1],s), "style":"fill:none;stroke:%s;stroke-width:%f;"%(color,width),"comment":str(comment)} )
		else :
			inkex.etree.SubElement( group, inkex.addNS('path','svg'), {"d": "m %s,%s l 10,10 -20,-20 10,10 -10,10, 20,-20"%(x[0],x[1]), "style":"fill:none;stroke:%s;stroke-width:%f;"%(color,width),"comment":str(comment)} )
	
			
										
################################################################################
###
###		Effect
###
###		Main function
###
################################################################################
	def effect(self) :
		if (self.options.create_bearing_gear and self.options.active_tab == '"help"' ):
			self.error(
"""English spport forum:
http://www.cnc-club.ru/forum/viewforum.php?f=33		
			
Russian support forum:
	http://cnc-club.ru/forum/viewtopic.php?f=15&t=287
			
			
			""")		
			return	
				
		time_ = time.time()
		self.distance = self.options.distance * self.options.units
		self.options.linear_distance *= self.options.units
		self.options.gear_r *= self.options.units
		self.options.bearings_d *= self.options.units
		
		self.first_teeth = self.options.selected_puley_teeth
		self.second_teeth = self.options.generated_puley_teeth
		self.num = self.options.number_of_copies

		def add_circle(c,r, center_mark=False) :
			return (
					"M%s,%s a %s,%s 0 1 1 %s,0 %s,%s 0 1 1 %s,0 z " % ( c[0]+r,c[1], r,r, -2*r, r,r, 2*r) + 
						#c+r,c   r r        -2*r  r   r       2*r
					( ("M %s,%s l -11,-11 22,22 -11,-11 -11,11, 22,-22" % (c[0],c[1])) if center_mark else "")
				)	
				
					
	
		c = list(self.view_center)
		d = ""
		d1 = ""		
		if (self.options.create_bearing_gear and self.options.active_tab == '"linear_shaft"' ):
			r = (self.options.gear_r*2+self.options.bearings_d/2+100)
			for i in range(self.options.number_of_bearings):
				a = i*math.pi*2/self.options.number_of_bearings
				d += add_circle([c[0]+math.cos(a)*self.options.gear_r, c[1]+math.sin(a)*self.options.gear_r	], self.options.bearings_d/2)
				d1 += add_circle([c[0]+math.cos(a)*self.options.gear_r, -r+c[1]+math.sin(a)*self.options.gear_r	], self.options.bearings_d/2, True)
			d1 += add_circle([c[0], c[1]-r], self.options.gear_r, True)


			path = inkex.etree.SubElement( self.current_layer, inkex.addNS('path','svg'), {"d":d, "style": "stroke:#4d4d4d;fill:#ececec"} )
			path1 = inkex.etree.SubElement( self.current_layer, inkex.addNS('path','svg'), {"d":d1, "style": "stroke:#4d4d4d;fill:#ececec"} )							

		
			csp = cubicsuperpath.parsePath(d)
			minx,miny,maxx,maxy = csp_true_bounds(csp)
			c1 = [ (maxx[0] + minx[0])/2, (maxy[1] + miny[1])/2 ] 
			path.set(inkex.addNS('transform-center-x','inkscape'),str(-c1[0]+c[0]))
			path.set(inkex.addNS('transform-center-y','inkscape'),str(-c1[1]+c[1]))


		
			self.selected = {"1":path}
		
		
		if len(self.selected)!=1: self.error('Please select one and only one path!','error')
		path = self.selected.values()[0]
		if "d" not in path.keys() :  self.error('Please select a "Path"!','error')
		csp = cubicsuperpath.parsePath(path.get("d"))
		
		center_group = inkex.etree.SubElement( path.getparent(),  inkex.addNS('g','svg') )
		group = inkex.etree.SubElement( path.getparent(),  inkex.addNS('g','svg') )
		attrib = path.attrib

		if "transform" in path.keys() :	
				t = path.get('transform')
				t = simpletransform.parseTransform(t)
				simpletransform.applyTransformToPath(t,csp)
				path.set("transform", "")
				path.set('d', cubicsuperpath.formatPath(csp))
				inkex.etree.SubElement( group, inkex.addNS('path','svg'), {"d": cubicsuperpath.formatPath(csp)})							
				self.error(path.get("transform"))

		# Will have to find center of bbox
		minx,miny,maxx,maxy = csp_true_bounds(csp)

		c1 = [ (maxx[0] + minx[0])/2, (maxy[1] + miny[1])/2 ] 
		w,h  = max(maxx[0]-minx[0], abs(c1[0]-minx[0]), abs(c1[0]-maxx[0]) ), max(maxy[1]-miny[1], abs(c1[1]-miny[1]), abs(c1[1]-maxy[1]) )
		if inkex.addNS('transform-center-x','inkscape') in path.keys() : c1[0] += float(path.get(inkex.addNS('transform-center-x','inkscape')))
		if inkex.addNS('transform-center-y','inkscape') in path.keys() : c1[1] -= float(path.get(inkex.addNS('transform-center-y','inkscape')))
		
		c2 = [c1[0]-self.distance, c1[1]]

		def get_transform(c1,c2,a1,a2):
			[c1x,c1y], [c2x,c2y] = c1,c2
			# move c1 to 0 
			transform = [ [1,0,-c1x], [0,1,-c1y] ]
			# Rotate to a1 to 0 
			transform = simpletransform.composeTransform([ [math.cos(a1), -math.sin(a1), 0], [math.sin(a1), math.cos(a1), 0] ], transform )
			# Move c2 to 0 			
			transform = simpletransform.composeTransform( [ [1,0,-c2x+c1x], [0,1,-c2y+c1y] ], transform)
			# Rotate to a2 to 0 
			transform = simpletransform.composeTransform( [ [math.cos(a2), -math.sin(a2), 0], [math.sin(a2), math.cos(a2), 0] ] , transform)
			# move c2 back to c2
			transform = simpletransform.composeTransform([ [1,0,c2x], [0,1,c2y] ], transform)
			return transform
		if not self.options.active_tab == '"linear_shaft"' :
			# Radial gear
			if not self.options.variable_speed :
				for i in range(self.num):
					
					alpha = math.pi*2*(i)/self.num*self.options.rev
					alpha1 = alpha*self.second_teeth
					alpha2 = alpha*self.first_teeth 

					transform = get_transform(c1,c2,alpha1,alpha2)
					# add path's clone 
					attrib["transform"] = simpletransform.formatTransform(transform)
					inkex.etree.SubElement( group, inkex.addNS('path','svg'), attrib )							

			else :
				def get_k (csp, c1,d):
					c2 = [c1[0]-d, c1[1]]
					alpha2_ = []
				
					for n in range(self.num):
						alpha1 = math.pi*2*(n)/self.num * self.second_teeth
						d_alpha1 = math.pi*2/self.num * self.second_teeth
						csp_ = [[[p[:] for p in point] for point in subpath] for subpath in csp ]
						transform = get_transform(c1,c2,alpha1,0)
						simpletransform.applyTransformToPath(transform, csp_)
						# r2 = distance to second gear's center
						[r2, i,j,t] = csp_to_point_distance(csp_, c2, dist_bounds = [0,1e100], tolerance=.001)
						r2 = math.sqrt(r2)

						p = csp_at_t(csp_[i][j-1],csp_[i][j],t)
						r1 = math.sqrt((p[0]-c1[0])**2 +(p[1]-c1[1])**2)
						# d_alpha2 = rotation speed factor
						if r2 == 0 : return 1e100, []
						alpha2_.append(d_alpha1 * r1/r2)
					return math.pi*2 * self.first_teeth / sum(alpha2_), alpha2_
				
				
				#get K, firs calculate aproximate turn and then get the K
				if not self.options.optimize_distance : 
					K, alpha2_ = get_k(csp,c1,self.distance)
				else : 		
					#set min and max distances
					dists = [0., math.sqrt(w**2+h**2)*(1+self.second_teeth/self.first_teeth) ]
					first = True
					for i in range(optimize_distance_max_iters):
						d = (dists[0]+dists[1])/2 #	if not first and not dists[0]<self.distance<dists[1] else self.distance
						first = False
						K, alpha2_ = get_k(csp,c1,d)
						if K>1 : 
							dists = [dists[0], d]
						else: 
							dists = [d, dists[1]] 
						if abs(1.-K)< optimize_distance_tolerance : 
							break
						#self.error(str((i,K,d)))	
					self.distance = d			
				c2 = [c1[0]-self.distance, c1[1]]	
				# Now create the paths
				alpha2 = 0
				for i in range(self.num):
					alpha = math.pi*2*(i)/self.num
					alpha1 = alpha*self.second_teeth
					alpha2 += alpha2_[i] * K

				
					transform = get_transform(c1,c2,alpha1,alpha2)
					# add path's clone 
					attrib["transform"] = simpletransform.formatTransform(transform)
					inkex.etree.SubElement( group, inkex.addNS('path','svg'), attrib )							
							
				self.error("Optimized distance: %s. "%(self.distance/self.options.units))
				self.error("Optimized parameter K: %s (the closer to 1.0 the better)."%K)
				self.error("Time elapsed %s seconds." %(time.time()-time_))			
			self.draw_pointer(c1, color = "#000080", width = 1, group = center_group )
			self.draw_pointer(c2, color = "#000080", width = 1, group = center_group )
		else :
			# Linear gear
			c1x,c1y = c1[0], c1[1]
			i = 0 
			self.distance = self.options.linear_distance
			if self.options.standatd_distance and self.options.create_bearing_gear: 
				self.distance = self.options.gear_r+self.options.bearings_d/2
				
			while i<=self.options.rev*self.num :
				a = math.pi*2*i/self.num
				d = self.distance * a 
				# Move c to 0  			
				transform = [ [1,0,-c1x], [0,1,-c1y] ]
				# Rotate to a
				transform = simpletransform.composeTransform( [ [math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0] ] , transform)
				# Move 0 to c + d  			
				transform = simpletransform.composeTransform( [ [1,0,c1x+d], [0,1,c1y] ], transform)
				attrib["transform"] = simpletransform.formatTransform(transform)
				inkex.etree.SubElement( group, inkex.addNS('path','svg'), attrib )							
				i += 1

e = Unordinary_gears()
e.affect()					

