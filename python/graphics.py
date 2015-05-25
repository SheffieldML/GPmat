from pylab import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys
import time
from Tkinter import *

#Creates a graphical representation of the data
class Graphics() :

    def __init__(self,dt,points,limits,learner,reader,animate=False):

        self.dt = dt
        self.points = points
        self.limits = limits
        self.drawn = False
        self.reader = reader
        self.learner = learner
        self.animate = animate
        #Create openGL window
        name = 'MoCap ToolBox'
        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(800,800)
        glutCreateWindow(name)

        glShadeModel(GL_SMOOTH)
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glutDisplayFunc(self.animateBegin)
        glMatrixMode(GL_PROJECTION)

        camerax,cameray,cameraz,zmiddle = self.findCamera()
        perspective = (cameraz-self.limits[5])+50

        gluPerspective(45,1.,1.,perspective)
        glMatrixMode(GL_MODELVIEW)

        #Lights
        self.setUpLights()

        #Camera

        gluLookAt(camerax,cameray,cameraz,
                  camerax,cameray,zmiddle,
                  0,1,0)

        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

        #Action!
        glutMainLoop()

    def setUpLights(self):

        lightPosition = [0,0]

        for i in range(2):
            lightPosition[i] =[self.limits[square((i/4)+1)],self.limits[((mod(i,4)>1)*3)+2],self.limits[mod(i,2)*3],0.]

        self.light0(lightPosition[0])
        self.light1(lightPosition[1])

    def findCamera(self):
        dist = zeros(2)
        dist[0] = self.limits[1]-self.limits[4]
        dist[1] = self.limits[0]-self.limits[3]

        cameraz = (dist.max())/(2*math.tan(pi/8))+self.limits[2]
        cameray = (self.limits[1] + self.limits[4])/2.
        camerax = (self.limits[0] + self.limits[3])/2.

        zmiddle = (self.limits[2]+self.limits[5])/2
        return camerax,cameray,cameraz,zmiddle


    def light0(self,position):
        lightZeroColor = [0.8,1.0,0.8,1.0] #green tinged
        glLightfv(GL_LIGHT0, GL_POSITION, position)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor)
        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1)
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05)
        glEnable(GL_LIGHT0)

    def light1(self,position):
        lightZeroColor = [0.8,1.0,0.8,1.0] #green tinged
        glLightfv(GL_LIGHT1, GL_POSITION, position)
        glLightfv(GL_LIGHT1, GL_AMBIENT, lightZeroColor)
        glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.1)
        glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.05)
        glEnable(GL_LIGHT1)

    def animateBegin(self):

        if(self.drawn==False):
            glMaterialfv(GL_FRONT,GL_DIFFUSE,[1.0,0.,0.,1.])
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
            if self.animate:
                for points in self.points:
                    self.draw(points)
            else:
                self.draw(self.points[0])
                cs = CanvasScreen(self.learner,self.reader)
                cs.canvasDisplay(self)
            self.drawn = True


    def draw(self,points):
        glMaterialfv(GL_FRONT,GL_DIFFUSE,[1.0,0.,0.,1.])
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
        for i in range(len(points[0])):
            glPushMatrix()
            glTranslated(points[0][i],points[1][i],points[2][i])
            glutSolidSphere(1,20,20)
            glPopMatrix()
        glutSwapBuffers()
        return

class CanvasScreen():

    def __init__(self,learner,reader):

        self.SCREENSIZE = 500.
        self.learner = learner
        self.reader = reader

    def canvasDisplay(self,graphics):
        self.dimensions = []
        self.graphics = graphics

        master = Tk()
        self.box = Canvas(master, width=self.SCREENSIZE, height=self.SCREENSIZE)
        plotline = self.learner.getPath(0,1)

        self.plot(plotline)

        self.box.bind('<Motion>', self.motion)

        self.box.bind('<Button-1>',self.left)
        mainloop()

    def plot(self,line):
        self.box.delete(ALL)
        #Convert given line dependent on display size
        line = self.actual2Display(line)
        line = line.T
        #Convert from 2D array containing xs and ys into tuples of (x,y) coords
        z = [(line[0][a],line[1][a]) for a in range(len(line[0]))]
        #Plot the line
        self.box.create_line(z)
        self.box.pack()

    def motion(self,event):
        #Get current screen position, convert to actual value
        position = self.display2Actual([event.x, event.y])

        #Concat existing dimensions with the current mouse position
        dimensions = array(self.dimensions + position)

        #Get channels from learner
        channels = self.learner.getPose(dimensions,range(len(dimensions)))
        #print dimensions
        #Convert to matrix points
        oldPoints = self.reader.bvh2xyz(channels)

        #Split into plottable points
        points = self.reader.split(oldPoints)

        #Draw new pose in openGL
        self.graphics.draw(points)

    def left(self,event):

        #Concat the x coordinate of the click onto stored dimensions list
        self.dimensions += [event.x]

        #Get the path for the next two components
        path = self.learner.getPath(len(self.dimensions),len(self.dimensions)+1)
        self.plot(path)

    def actual2Display(self,x):
        self.min0 = min(x[:,0])*1.
        self.min1 = min(x[:,1])*1.

        print "Min X0"
        print self.min0
        print "Min X1"
        print self.min1

        self.max0 = max(x[:,0])*1.
        self.max1 = max(x[:,1])*1.

        print "Max X0"
        print self.max0
        print "Max X1"
        print self.max1

        #Set the minimum to 0
        x[:,0] -= self.min0
        x[:,1] -= self.min1

        self.span0 = max(x[:,0])*1.
        self.span1 = max(x[:,1])*1.

        #Set the max to screensize
        x[:,0] *= self.SCREENSIZE/self.span0
        x[:,1] *= self.SCREENSIZE/self.span1
        return x


    def display2Actual(self,x):

        b = [0,0]

        #NONLINEAR
        #b[0] = x[0]/30.
        #b[1] = x[1]/30.
        #LINEAR
        b[0] = (x[0]/(self.SCREENSIZE/(self.span0)))+self.min0
        b[1] = (x[1]/(self.SCREENSIZE/self.span1))+self.min1
        return b



