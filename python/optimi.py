import pdb
import numpy as np
import math

def create(trans):
    if trans=='exp':
        return expTransform()
    else:
        raise Exception('Unknown transform type.')
    

class optimisable:
    """Base class for optimisable models, implements gradient checks and various optimizers."""
    def __init__(self):
        if self.__class__ is optimisable:
            raise NotImplementedError, "optimisable is a base class and should not be created directly, only inherited"

    def getOptNumParams(self):
        pass
    def getOptParams(self):
        pass
    def setOptParams(self):
        pass

    def setVerbosity(self, val=2):
        self.verbosity = val
    def getVerbosity(self):
        return self.verbosity

    def gradients(self):
        pass
    def objective(self):
        pass
    def getParam(self):
        pass
    def setParam(self, parameters):
        pass

    def checkGradients(self):
        change = 1e-6
        param = self.getParam()
        origParam = param.copy()
        analyticGrad = self.gradient(param)
        numericalDiff = np.zeros(analyticGrad.shape)
        for j in range(len(param)):
            origVal = origParam[j]
            param[j] = origVal + change
            #self.setParam(param)
            objectivePlus = self.objective(param)
            param[j] = origVal - change
            #self.setParam(param)
            objectiveMinus = self.objective(param)
            numericalDiff[j] = 0.5*(objectivePlus - objectiveMinus)/change
            param[j] = origVal
        self.setParam(origParam)
        print "Numerical differences:, ", numericalDiff
        print "Analytic gradients: ", analyticGrad 
        print "Differences: ", numericalDiff - analyticGrad
        
        

class transformable(optimisable):
    """This class allows the parameters of a model to be transformed."""

    def __init__(self):
	optimisable.__init__(self)
        self.transArray = []
    
    def addTransform(self, trans, index):
        if isinstance(trans, basestring):
            trans = create(trans)
        self.transArray.append(transforms(trans, index))

    def getNumTransforms(self):
        return len(self.transArray)

    def getTransform(self, i):
        return self.transArray[i].trans

    def getTransformIndex(self, i):
        return self.transArray[i].index

    def getParam(self):
        return self.extractTransParam()

    def setParam(self, params):
        self.expandTransParam(params)

    def expandParam(self, params):
        pass

    def extractParam(self):
        pass

    def expandTransParam(self, params):
        newParams = params.copy()
        for i in range(self.getNumTransforms()):
            ind = self.getTransformIndex(i)
            newParams[ind] = self.getTransform(i).atox(newParams[ind])
        self.expandParam(newParams)

    def extractTransParam(self):
        params = self.extractParam()
        for i in range(self.getNumTransforms()):
            ind = self.getTransformIndex(i)
            params[ind] = self.getTransform(i).xtoa(params[ind])
        return params

    def gradTransParam(self):
        params = self.extractParam()
        g = self.gradParam(self)
        for i in range(self.getNumTransforms()):
            ind = getTransformIndex(i)
            facts = self.getTransform(i).gradfact(params[ind])
            g[ind] = g[ind]*facts
        return g

    def gradParam(self):
        pass


class transform:
    """Base class for the transforms."""

    def __init__(self):
        pass

    def atox(self, x):
        pass

    def xtoa(self, x):
        pass

    def gradfact(self, x):
        pass

class expTransform(transform):
    """EXPTRANSFORM Constrains a parameter to be positive through exponentiation.
    FORMAT
    DESC contains commands to constrain parameters to be positive via
    exponentiation.
    ARG x : input argument.
    ARG y : return argument.
    ARG transform : type of transform, 'atox' maps a value into
    the transformed space (i.e. makes it positive). 'xtoa' maps the
    parameter back from transformed space to the original
    space. 'gradfact' gives the factor needed to correct gradients
    with respect to the transformed parameter.
    
    SEEALSO : negLogLogitTransform, sigmoidTransform
    
    COPYRIGHT : Neil D. Lawrence, 2004-2007, 2009

    """
    
    def __init__(self):
        transform.__init__(self)

    def atox(self, a):
        limVal = 36
        x = np.zeros(a.shape)
        index = np.nonzero(a<-limVal)[0]
        x[index] = np.exp(-limVal)
        a[index] = np.nan
        index = np.nonzero(a<limVal)
        x[index] = np.exp(a[index])
        a[index] = np.nan
        index = np.nonzero(np.logical_not(np.isnan(a)))
        if len(index[0])>0:
            x[index] = np.exp(limVal)
        return x

    def xtoa(self, x):
        return np.log(x)
    
    def gradfact(self, x):
        y = x.copy()
        return y


def defaultConstraint(constraint):

    '''% OPTIMIDEFAULTCONSTRAINT Returns function for parameter constraint.
    % FORMAT
    % DESC returns the current default function for constraining a
    % parameter. Formerly (up to version KERN 0.163) this was
    % 'negLogLogit' for positive constraints, as this keeps things roughly linear in the
    % positive half space, however, it is more standard to use 'exp'
    % (i.e. optimise in the log space). This function allows you to
    % change the option globally. The defaults are given below.
    % ARG constraint : the type of constraint you want to place on the
    % parameter, options include 'positive' (gives an 'exp' constraint)
    % and 'zeroone' (gives a 'sigmoid' constraint).
    % RETURN str : the type of function used to apply the constraint
    % from the 'optimi' toolbox.
    %
    % SEEALSO : expTransform, sigmoidTransform, linearTransform, negLogLogitTransform
    % COPYRIGHT : Neil D. Lawrence, 2006, 2009
    
    % OPTIMI

    '''
    
    if constraint == 'positive':
        str = 'exp'
    elif constraint == 'zeroone':
        str = 'sigmoid'
    else:
        raise Exception("Unrecognized requested default constraint.")
    return str

class transforms:
    """Simple transforms class for containing the index of parameters to be transformed as well as the transform types."""

    def __init__(self, trans, index):
        self.index = index
        self.trans = trans

class fixable(transformable):
	"""A class to allow objects fo also have fixed values, as well as transformed ones"""
	def __init__(self):
		transformable.__init__(self)
		self.fixed_index = np.array([],dtype=np.int)
		self.fixed_values = np.array([],dtype=np.float64)

	def addFix(self,index,value):
		self.fixed_index = np.hstack((self.fixed_index,index))
		self.fixed_values = np.hstack((self.fixed_values,value))
		
		#set the model's value to the fixed value!
		p = self.extractTransParam()
		p[index] = value
		self.expandTransParam(p)

	def setParam(self,params):
		allparam = np.zeros(params.size+self.fixed_values.size)
		allparam[self.fixed_index] = self.fixed_values
		unfixed_index = np.nonzero(allparam==0)[0]
		allparam[unfixed_index] = params
		self.expandTransParam(allparam)

	def getParam(self):
		param = self.extractTransParam()
		unfixed_index = np.setdiff1d(np.arange(param.size),self.fixed_index)
		return param[unfixed_index]

class tieable(fixable):
	"""allows some of the parameters to be tied together"""
	def __init__(self):
		fixable.__init__(self)
		self.tied_visible_index = np.array([],dtype=np.int)
		self.tied_hidden_index = np.array([],dtype=np.int)
	
	def tieParam(self,index1,index2):
		self.tied_visible_index = np.hstack((self.tied_visible_index,index1))
		self.tied_hidden_index = np.hstack((self.tied_hidden_index,index2))
		#set tied parameters to their values... not sure if this python hack is healthy?
		p = fixable.getParam(self)
		p[index2] = p[index1]
		fixable.setParam(self,p)

	def setParam(self,params):
		old_param = fixable.getParam(self)
		N = old_param.size
		newparam = np.zeros(N)
		count = 0
		for i in range(N):
			if not i in self.tied_hidden_index:
				newparam[i]=params[count]
				count += 1
		newparam[self.tied_hidden_index] = newparam[self.tied_visible_index]
		pdb.set_trace()
		fixable.setParam(self,newparam)

	def getParam(self):
		p = fixable.getParam(self)
		return np.delete(p,self.tied_hidden_index,0)



	
