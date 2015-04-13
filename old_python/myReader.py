import cgkit.bvh as bvh
#from cgkit import *

from pylab import *
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class myReader(bvh.BVHReader) :
	def __init__(self, filename):
		self.screen = None
		self.filename = filename
		# A list of unprocessed tokens (strings)
		self.tokenlist = []
		# The current line number
		self.linenr = 0
		self.frame = 0
		# Root node
		self._root = None
		self._nodestack = []
		self.points = []
		# Total number of channels
		self._numchannels = 0
		self.count = 0
		self.nodes = []
		self.frames = 0
	
	def onMotion(self,frames,dt):
		self.screen = graphics(dt)
		self.frames = frames

	def onFrame(self,channels):		
		self.points.append(self.bvh2xyz(channels))
		self.frame += 1
		if(self.frame ==self.frames):
			self.screen.limits = self.split()
			self.screen.setLimits()
			self.screen.draw(self.points)
			
	def recurse(self,node):
		if(node.children) :
			for n in node.children:
				self.nodes.append(n)				
				self.count+=1
				self.recurse(n)

	def onHierarchy(self,root):
		self.count+=1
		self.nodes.append(root)
		self.recurse(root)
		self.results = [[0,0]]
		for n in range(len(self.nodes)):
			for x in self.nodes[n].children:	
				child = self.findNode(x)
				self.results.append([child,n])
		self.mySort()	
				
			
	def findNode(self,node):
		for n in range(len(self.nodes)):
			if(self.nodes[n] == node):
				return n
	
	def mySort(self):
		temp = []
		
		for n in self.results :
			if len(temp)==0 :
				temp.append(n)
			else :
				for i in range(len(temp)) :
					if temp[i][0]> n[0]   :
						temp.insert(i,n)
						break
				else:
					temp.append(n)
		self.results = temp				

	def bvh2xyz(self,channels):
		curChan = 0
		xyzStruct = []
		
		for i in range(len(self.nodes)):
			if (len(self.nodes[i].channels) == 6) :
				xpos = channels[curChan]
				ypos = channels[curChan+1]
				zpos = channels[curChan+2]
				curChan+=3
			else:
				xpos = 0
				ypos = 0
				zpos = 0
			if(len(self.nodes[i].channels) >2):
				zangle = math.radians(channels[curChan])
				xangle = math.radians(channels[curChan+1])
				yangle = math.radians(channels[curChan+2])
				curChan+=3
			else :
				xangle = 0.
				yangle = 0.
				zangle = 0.
			offsets = matrix('0.,0.,0.')
			thisRotation = self.rotationMatrix(xangle,yangle,zangle)
			thisPosition = [xpos, ypos, zpos]
			thisPosition = offsets + thisPosition
			struct = position()
			xyzStruct.append(struct)
			offsets += self.nodes[i].offset
			if i==0 :
				xyzStruct[i].position = offsets + thisPosition
				xyzStruct[i].rotation = thisRotation
			else :
				parent = self.results[i][1]
				xyzStruct[i].position = (offsets + thisPosition)*xyzStruct[parent].rotation + xyzStruct[parent].position
				xyzStruct[i].rotation = thisRotation*xyzStruct[parent].rotation
		points = []		
		
		for m in xyzStruct :
			points.append(m.position)
		return points

	def rotationMatrix(self,xangle,yangle,zangle,order='zxy'):
		c1 = math.cos(xangle) 
		
		c2 = math.cos(yangle) 
		c3 = math.cos(zangle) 
		s1 = math.sin(xangle)
		s2 = math.sin(yangle)
		s3 = math.sin(zangle)
		rM = array([[c2*c3-s1*s2*s3, c2*s3+s1*s2*c3, -s2*c1],
            		 [-c1*s3, c1*c3, s1],
            		 [s2*c3+c2*s1*s3, s2*s3-c2*s1*c3, c2*c1]])
		return matrix(rM)
			
	def read(self):
		bvh.BVHReader.read(self)
		#return self.points
		
	def split(self):
		length = len(self.points[0])
		size = len(self.points)
		newPoints = zeros([size,3,length],'float')
		zmax = xmax = ymax = 0
		zmin = xmin = ymin = 0
		for i in range(size):
			zs = zeros(length)
			xs = zeros(length)
			ys = zeros(length)
			for j in range(0,length):
				
				line = self.points[i][j].getA()
				zs[j] = line[0][0]
				xs[j] = line[0][1]
				ys[j] = line[0][2]
				
			if(i==0):
				zmax = max(zs)
				xmax = max(xs)
				ymax = max(ys)
			maxz = max(zs)
			maxy = max(ys)
			maxx = max(xs)
			if(maxz > zmax):
				zmax = maxz
			if(maxx > xmax):
				xmax = maxx
			if(maxy > ymax):
				ymax = maxy
				
			
			minz = min(zs)
			miny = min(ys)
			minx = min(xs)
			if(minz < zmin):
				zmin = minz
			if(minx < xmin):
				xmin = minx
			if(miny < ymin):
				ymin = miny
			newPoints[i][0] = zs
			newPoints[i][1] = xs
			newPoints[i][2] = ys
		
		self.points = newPoints
		return ([zmax,xmax,ymax,zmin,xmin,ymin])
			
class position():
	def __int__(self,rotation,position):
		self.rotation = rotation
		self.position = position

class graphics() :
	def __init__(self,dt):
		self.dt = dt
		self.fig = plt.figure()
		self.ax = Axes3D(self.fig)
		#self.fig.add_subplot(111, projection='3d')
		self.ax.set_xlabel('X Label z axis')
		self.ax.set_ylabel('Y Label x axis')
		self.ax.set_zlabel('Z Label Y Axis')
		self.limits = zeros(6)
		#self.ax.autoscale(enable=False,axis='both',tight=False)
	
	def setLimits(self):		
		self.ax.set_ylim3d(self.limits[3],self.limits[0])
		self.ax.set_zlim3d(self.limits[4],self.limits[1])
		self.ax.set_xlim3d(self.limits[5],self.limits[2])
		
		plt.show()
		
	def draw(self,data):
		print data[0]
		for points in data :
			
			#print "draw"
			#o = time.time()
			self.ax.set_autoscale_on(False)
			self.ax.scatter3D(points[2],points[0],points[1],c='r',marker='o')
			draw()
			#f = time.time()
			
		#	time.sleep(self.dt)
			self.ax.clear()
			#print f-o






















