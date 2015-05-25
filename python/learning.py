from pylab import *
import copy
import GPy

class abstractLearner():

    def __init__(self,fullChannelList):
        self.channels = fullChannelList

    def getPose(self,dimensions):
        #Needs to return a pose of an equal number of channels
        #From a single frame of fullChannelList
        pass

    def getPath(component1,component2):
        #components are indices of components
        #return a list of 2D coordinates in the given components
        #from the reduced feature space
        pass

    def getPoint(self,index,component1,component2):
        # Gets the reduced point of a particular index
        pass

class linearLearner(abstractLearner):

    def __init__(self,fullChannelList):
        #Import training data
        self.channels = array(fullChannelList)
        #Compute covariance
        s = cov(self.channels.T)
        #Get eigenvalues
        self.u,self.v = eig(s)


    def getPath(self,component1,component2):
        #Project the sequence onto two prinicpal components
        #In order to visualise
        w = self.getW([component1,component2])
        #return the path for the given components
        x = array(dot(self.channels,w))

        return x

    def getPose(self,dimensions,componentIndices):
        w = self.getW(componentIndices)
        newChannels = dot(w,dimensions) + self.channels.mean(axis=0)
        return newChannels

    def getW(self,componentIndices):
        #make w from relevant components
        w = []
        for i in componentIndices:
            w += [real(self.v[:,i])]
        return array(w).T

    def getPoint(self,index, component1,component2):
        d = self.getPath(component1,component2)
        return d[index]

# This line would produce the same projection you already have. But we will need a line like this when we do a basis function model.
        #w = dot(dot(inv(dot(x.T,x)),x.T),y)
        #print array(x).T[0]

class nonLinearLearner(abstractLearner):

    def __init__(self,fullChannelList):
        self.channels = array(fullChannelList)

        XPCA,W = GPy.util.linalg.PCA(self.channels,2)
        k = GPy.kern.rbf(XPCA, alpha=10.,gamma=1) + GPy.kern.white(XPCA,alpha=1.)
        k.constrain_positive('cmp')
        self.m = GPy.models.GPLVM(self.channels - self.channels.mean(axis=0),k)
        #self.m = GPy.models.simple_GP(self.channels - self.channels.mean(axis=0),XPCA,k)
        self.m.optimize(max_f_eval=1)

    def getPath(self,component1,component2):
        return np.copy(self.m.X)

    def getPose(self,dimensions,componentIndices):
        newDim = array(dimensions)

        a = self.m.predict(array([newDim]))
        return (a + self.channels.mean(axis=0)).T

    def getPoint(self,index,_1,_2):
        return np.copy(self.m.X[index])

















