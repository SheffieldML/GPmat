import math
import numpy as np

##SETUP
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'optimi', 'python'))
#sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'mlopy', 'netlab'))
import pdb
##ENDSETUP

#import netlab
import optimi
from sim import simComputeH, lnDiffErfs

def dist2(x1,x2):
	return np.sum(np.square(x1[:,np.newaxis,:]-x2[np.newaxis,:,:]),-1)

def test(kernType, numIn=4, tieParamNames=None):
    """% TEST Run some tests on the specified kernel.
    % FORMAT
    % DESC runs some tests on a kernel with the specified type to ensure it is
    % correctly implemented.
    % ARG kernType : type of kernel to test. For example, 'rbf' or
    % {'cmpnd', 'rbf', 'lin', 'white'}.
    % ARG numIn : the number of input dimensions (default is 4).
    % ARG tieParamNames : cell array of regular expressions for parameter
    % names that should be tied (default is none).
    % RETURN kern : the kernel that was generated for the tests.
    % 
    % FORMAT
    % DESC runs some tests on a given kernel structure to ensure it is
    % correctly implemented.
    % ARG kern : kernel structure to test.
    % ARG numIn : the number of input dimensions (default is 4).
    % ARG tieParamNames : cell array of regular expressions for parameter
    % names that should be tied (default is none).
    % RETURN kern : the kernel as it was used in the tests.
    % 
    % SEEALSO : create
    %
    % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2007, 2009
    %
    % MODIFICATIONS : Antti Honkela, 2007
    %
    % MODIFICATIONS : David Luengo, 2009
    
    """

    numData = 20

    if isinstance(kernType, kern):
        ker = kernType
        kernType = ker.type
        params = ker.extractTransParam()
        paramnames = ker.extractParamNames()
        paramExpand = np.eye(len(paramnames))
        paramPack = paramExpand
        toRemove = []
        kerCreate = False
    else:
        ker = create(kernType, numIn)
        kerCreate = True

    if tieParamNames is None:
        tieParamNames = []


    # For convolutional kernels starting at t=0 it does not make sense
    # to use negative inputs...
    if ker.positiveTime:
        x = np.abs(np.random.randn(numData, numIn))
        x2 = np.abs(np.random.randn(floor(numData/2), numIn))
    else:
        x = np.random.randn(numData, numIn)
        x2 = np.random.randn(numData, numIn)

    if kerCreate:
        if isinstance(ker, component): 
            for i in range(len(ker.comp)):
                if float(np.random.rand(1))>0.5:
                    indices = np.random.permutation(numIn)
                    indices = indices[0:int(np.ceil(np.random.rand(1)*numIn))]
                    ker.setIndex(i, indices)

        # Set the parameters randomly.
        params = ker.extractTransParam()
        paramNames = ker.extractParamNames()
        if not isinstance(tieParamNames, np.ndarray):
            tieParams = []
            paramExpand = np.eye(len(params))
            toRemove = []
            for l in range(len(tieParamNames)):
                #ties = strfind(paramnames, tieParamNames[l];
                ties = paramNameRegularExpressionLookup(model, tieParamNames[0])
                ties = regexp(paramnames, tieParamNames[l]);
                tieParams[0] = []
                for k in range(length(ties)):
                    if len(ties[k])>0:
                        tieParams[l].append(k)
                if len(tieParams)>0:
                    paramExpand[:, tieParams[0][0]] = sum(paramExpand[:, tieParams[0]], 1)
                    toRemove = np.array([toRemove, tieParams[l][1:end]])
            oldParamExpand = paramExpand
            paramExpand = np.zeros((paramExpand.shape[0], 
                                    paramExpand.shape[1]-len(toRemove)))
            count = 0
            for i in range(paramExpand.shape[1]):
                if not i in toRemove:
                    paramExpand[:, count] = oldParamExpand[:, i]
                    count += 1
        else:
            paramExpand = tieParamNames
            toRemove = []
            for k in range(paramExpand.shape[1]):
                I = np.nonzero(paramExpand[:, k])
                toRemove.append(I[1:])
        paramPack = paramExpand.T / sum(paramExpand.T, 1)
        params = np.dot(params, paramPack.T)
        params = np.random.standard_normal(params.shape)/np.abs(np.random.standard_normal(params.shape))
        ker.expandTransParam(np.dot(params, paramExpand.T))

    # Test for positive definiteness
    K = ker.compute(x)[0]
    e = np.linalg.eigvals(K)
    if min(e) > 0:
        print 'The kernel is positive definite.'
    else:
        print 'The kernel is not positive definite: max eig ', np.max(e), \
            ', min eig ', np.min(e)

    covGrad = np.ones(K.shape)
    epsilon = 1e-6
    Lplus = np.zeros(params.shape)
    Lminus = np.zeros(params.shape)
    params = np.dot(ker.extractTransParam(), paramPack.T)
    origParams = params.copy()
    for i in range(len(params)):
        params[i] = origParams[i] + epsilon
        ker.expandTransParam(np.dot(params,paramExpand.T))
        Lplus[i] = ker.compute(x)[0].sum()
        params[i] = origParams[i] - epsilon
        ker.expandTransParam(np.dot(params, paramExpand.T))
        Lminus[i] = ker.compute(x)[0].sum()
        params[i] = origParams[i]

    ker.expandTransParam(np.dot(params, paramExpand.T))
    names = ker.extractParamNames()
    toRemove.sort(reverse=True)
    for i in toRemove:
        del names[i]
    gLDiff = .5*(Lplus - Lminus)/epsilon
    g = np.dot(ker.gradTransParam(x, covGrad=covGrad), paramExpand)
        
    paramMaxDiff = np.abs(gLDiff-g).max()
    if paramMaxDiff > 2*epsilon:
        l = 0
        for i in range(len(names)):
            if l < len(names[i]):
                l = len(names[i])
  
        print ' '*l, '\tanalytic   diffs     delta'
        for i in range(len(names)):
            spaceLen = l - len(names[i])
            space = ' '*spaceLen
            print space, names[i], ':\t%4.6f\t%4.6f\t%4.6f' % \
                (g[i], gLDiff[i], gLDiff[i] - g[i])
#    try:
    Lplus = np.zeros(x.shape)
    LplusDiag = np.zeros(x.shape)
    Lminus = np.zeros(x.shape)
    LminusDiag = np.zeros(x.shape)
    gx = np.zeros(x.shape)
    gxDiag = np.zeros(x.shape)
    origX = x.copy()
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            x[i, j] = origX[i, j] + epsilon
            K = ker.compute(x)[0]
            Lplus[i, j] =  K.sum()
            LplusDiag[i, j] = K.trace()
            x[i, j] = origX[i, j] - epsilon
            K = ker.compute(x)[0]
            Lminus[i, j] = K.sum()
            LminusDiag[i, j] = K.trace()
            x[i, j] = origX[i, j]
        gx[i:i+1, :] = 2*np.sum(ker.gradX(x[i:i+1, :], x)[:, :, 0], 0)
        gxDiag[i:i+1, :] = ker.diagGradX(x[i:i+1, :])

    gXDiff = 0.5*(Lplus - Lminus)/epsilon
    xMaxDiff = abs(gx-gXDiff).max()

    if xMaxDiff > 2*epsilon:
        print 'gX'
        print gx
        print 'gXDiff'
        print gXDiff

    gXDiagDiff = .5*(LplusDiag - LminusDiag)/epsilon
    xDiagMaxDiff = abs(gxDiag-gXDiagDiff).sum()

    if xDiagMaxDiff > 2*epsilon:
        print 'gxDiag'
        print gxDiag
        print 'gXDiagDiff'
        print gXDiagDiff

 #   except Exception:
 #       print 'kern.gradX has an error.'
 #     #warning(lasterr)
 #       xMaxDiff = 0
 #       xDiagMaxDiff = 0

    K = ker.compute(x)[0]
    traceK =  K.trace()
    K2 = ker.diagCompute(x)
    traceK2 = K2.sum()
    traceDiff = traceK - traceK2
    if abs(traceDiff) > 2*epsilon:
        print 'kernDiagCompute is not in sync with kernCompute.'
        print 'diag(kernCompute)\tkernDiagCompute'
        print np.array([np.diag(K), K2])

    covGrad = np.ones((x.shape[0], x2.shape[0]))
    params = np.dot(ker.extractTransParam(),paramPack.T)
    origParams = params.copy()
    Lplus = np.zeros(params.shape)
    Lminus = np.zeros(params.shape)
    for i in range(len(params)):
        params[i] = origParams[i] + epsilon
        ker.expandTransParam(np.dot(params, paramExpand.T))
        Lplus[i] = ker.compute(x, x2)[0].sum()
        params[i] = origParams[i] - epsilon
        ker.expandTransParam(np.dot(params, paramExpand.T))
        Lminus[i] = ker.compute(x, x2)[0].sum()
        params[i] = origParams[i]

    ker.expandTransParam(np.dot(params, paramExpand.T))
    names = ker.extractParamNames()
    toRemove.sort(reverse=True)
    for i in toRemove:
        del names[i]
    gL2Diff = .5*(Lplus - Lminus)/epsilon
    g = np.dot(ker.gradTransParam(x, x2, covGrad=covGrad), paramExpand)
    
    param2MaxDiff = abs(gL2Diff-g).max()
    if param2MaxDiff > 2*epsilon:
        l = 0
        for i in range(len(names)):
            if l < len(names[i]):
                l = len(names[i])

        print ' '*l, '\tanalytic   diffs     delta'
        for i in range(len(names)):
            spaceLen = l - len(names[i])
            space = ' '*spaceLen
            print space, names[i], ':\t%4.6f\t%4.6f\t%4.6f' % \
                (g[i], gL2Diff[i], gL2Diff[i] - g[i])
    #pause(0)

    print 'Trace max diff: %2.6f.' % (traceDiff)
    print 'Param max diff: %2.6f.' % (paramMaxDiff)
    print 'Param X2 max diff: %2.6f.' % (param2MaxDiff)
    print 'X max diff: %2.6f.' % (xMaxDiff)
    print 'XDiag max diff: %2.6f.' % (xDiagMaxDiff)
    
    ker.display()
    return ker

    #% We don't test kernCompute(kern, x, x2) here at all!



def create(type, inDim=None, X=None):
    """CREATE Helper function for creating a kernel out of a list of kernel
    types."""
    if isinstance(type, basestring):
        if type == 'rbf':
            return rbf(inDim, X)
        elif type == 'lin':
            return lin(inDim, X)
        elif type == 'bias':
            return bias(inDim, X)
        elif type == 'white':
            return white(inDim, X)
        else:
            raise Exception('Unknown kernel type in create.')
    else:
        if type[0] == 'cmpnd':
            ker = cmpnd(inDim, X)

        #elif type[0] == 'tensor':
        #    ker = tensor(inDim, X)
        
        else:
            ker = cmpnd(inDim, X)
            ker.addKern(create(type[0], inDim, X))

        for i in range(1, len(type)):
            ker.addKern(create(type[i], inDim, X))

        return ker


class kern(optimi.tieable):
    """The base kernel class from which all kernels are derived."""

    def __init__(self, inDim=None, X=None):
        optimi.tieable.__init__(self)
        self.whiteVariance = None # whether the covariance has a white
                                  # noise component
        self.positiveTime = False # default value, set true for
                                  # laplace convolutional kernels
                                  # starting from zero.
        if inDim is None and X is not None:
            inDim = X.shape[1]
        elif inDim is None:
            # if input dimension not specified assume it is 1.
            inDim = 1
        self.inputDimension = inDim
        self.paramInit(inDim)

    def paramInit(self, inDim=None):
        pass

    def compute(self, x, x2=None):
        pass

    def diagCompute(self, x):
        pass

    def diagGradX(self, X):
        pass

    def diagGradient(self, X, covDiag):
        g = np.zeros(self.nParams)
        for i in range(X.shape[0]):
            g += self.gradient(x[i], covGrad=covDiag[i])
        return g

    def display(self, numSpaces=0):
        pass

    def gradX(self, X, X2):
        pass

    def gradient(self, X, X2=None, covGrad=None):
        pass

    def gradTransParam(self, X, X2=None, covGrad=None):
        """Return the gradient of the transformed parameters."""
        params = self.extractParam()
        g = self.gradient(X, X2=X2, covGrad=covGrad)
        for i in range(self.getNumTransforms()):
            ind = self.getTransformIndex(i)
            facts = self.getTransform(i).gradfact(params[ind])
            g[ind] = g[ind]*facts
        return g


class component(kern):
    """The base class for kernels which have several components to
    them (such as compound or tensor kernels).

    """

    def __init__(self, inDim=None, X=None):
        kern.__init__(self, inDim, X)
        self.paramGroups = np.array([])
        self.comp = []
        self.index = []
        self.nParams = 0

    def addKern(self, ker):
        oldNParams = self.nParams
        self.comp.append(ker)
        self.index.append(None)
        self.nParams += ker.nParams
        for i in range(ker.getNumTransforms()):
            ind = ker.getTransformIndex(i)
            for j in range(len(ind)):
                ind[j] += oldNParams
            self.addTransform(ker.getTransform(i), ind)      
            self.stationary = self.stationary and ker.stationary
        self.paramGroups = np.eye(self.nParams)
        return len(self.comp)-1

    def setIndex(self, component, indices):
        """ SETINDEX Set the indices in the compound kernel.
        
        """

        if indices.ndim > 1:
            raise Exception('Indices should be one dimensional.')
        if np.max(indices) > self.inputDimension:
            raise Exception('Indices are larger than kernel input dimension')
        elif np.min(indices) < 0:
            raise Exception('Indices should be greater than zero.');
  
        self.comp[component].inputDimension = len(indices);
        self.index[component] = indices
        self.comp[component].paramInit(self.comp[component])
        for i in range(len(self.comp)):
            self.nParams = self.nParams + self.comp[i].nParams;


        self.paramGroups = np.eye(self.nParams)


class cmpnd(component):
    """The is short for compound kernel and it is developed
    by summing several kernels together."""
    def __init__(self, inDim=None, X=None):
        component.__init__(self, inDim, X)

    def paramInit(self, inDim=None, X=None):
        self.type = 'cmpnd'
        self.stationary = True

    def compute(self, x, x2=None):
        if x2 is not None:
            i = 0
            if self.index[i] is not None:
                # only part of the data is involved in the kernel.
                k = self.comp[i].compute(x[:, self.index[i]], 
                                         x2[:, self.index[i]])[0]
            else:
                # all the data is involved with the kernel.
                k = self.comp[i].compute(x, x2)[0]
            for i in range(1, len(self.comp)):
                if self.index[i] is not None:
                    # only part of the data is involved in the kernel.
                    k += self.comp[i].compute(x[:, self.index[i]], 
                                              x2[:, self.index[i]])[0]
                else:
                    # all the data is involved with the kernel.
                    k += self.comp[i].compute(x, x2)[0]
        else:
            i = 0
            if self.index[i] is not None:
                # only part of the data is involved with the kernel.
                k = self.comp[i].compute(x[:, self.index[i]])[0]
            else:
                # all the data is involved with the kernel.
                k = self.comp[i].compute(x)[0]
            for i in range(1, len(self.comp)):
                if self.index[i] is not None:
                    # only part of the data is involved with the kernel.
                    k += self.comp[i].compute(x[:, self.index[i]])[0]
                else:
                    # all the data is involved with the kernel.
                    k += self.comp[i].compute(x)[0]
        return (k, )

    def diagCompute(self, x):
        
        """% DIAGCOMPUTE Compute diagonal of CMPND kernel.
        % FORMAT
        % DESC computes the diagonal of the kernel matrix for the
        % compound kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : a vector containing the diagonal of the kernel matrix
        % computed at the given points.
        %
        % SEEALSO : cmpndKernParamInit, kernDiagCompute, create, cmpndKernCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        i = 0
        if self.index[i] is not None:
            # only part of the data is involved with the kernel.
            k  = self.comp[i].diagCompute(x[:, self.index[i]])
        else:
            # all the data is involved with the kernel.
            k  = self.comp[i].diagCompute(x)
        for i in range(1, len(self.comp)):
            if self.index[i] is not None:
                # only part of the data is involved with the kernel.
                k += self.comp[i].diagCompute(x[:, self.index[i]])
            else:
                # all the data is involved with the kernel.
                k += self.comp[i].diagCompute(x)
        return k

    def diagGradX(self, X):        
        """Gradient of CMPND kernel's diagonal with respect to X.

        Computes the gradient of the diagonal of the compound kernel matrix with
        respect to the elements of the design matrix given in X.
        
        Arguments:
            X : the input data in the form of a design matrix.
        
        Return:
            gX : the gradients of the diagonal with respect to each element
            of X. The returned matrix has the same dimensions as X.
        
        See Also: 
        
        cmpndKernParamInit, kernDiagGradX, cmpndkernGradX
        
        Copyright : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """
        

        i = 0
        if self.index[i] is not None:
            # only part of the data is involved with the kernel.
            gX = np.zeros(X.shape)
            gX[:, self.index[i]] = self.comp[i].diagGradX(X[:, self.index[i]])
        else:
            # all the data is involved with the kernel.
            gX = self.comp[i].diagGradX(X)
        for i in range(1, len(self.comp)):
            if self.index[i] is not None:
                # only part of the data is involved with the kernel.
                gX[:, self.index[i]] = gX[:, self.index[i]] + self.comp[i].diagGradX(X[:, self.index[i]])
            else:
                # all the data is involved with the kernel.
                gX += self.comp[i].diagGradX(X)

        return gX

    def diagGradient(self, x, covDiag):
        """% DIAGGRADIENT Compute the gradient of the CMPND kernel's diagonal wrt parameters.
        % FORMAT
        % DESC computes the gradient of functions of the diagonal of the
        % compound kernel matrix with respect to the parameters of the kernel. The
        % parameters' gradients are returned in the order given by the
        % cmpndKernExtractParam command.
        % ARG kern : the kernel structure for which the gradients are
        % computed.
        % ARG x : the input data for which the gradient is being computed.
        % ARG factors : partial derivatives of the function of interest with
        % respect to the diagonal elements of the kernel.
        % RETURN g : gradients of the relevant function with respect to each
        % of the parameters. Ordering should match the ordering given in
        % cmpndKernExtractParam.
        %
        % SEEALSO : cmpndKernParamInit, kernDiagGradient, cmpndKernExtractParam, cmpndKernGradient
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """


        g = np.zeros(self.nParams)
        startVal = 0;
        endVal = 0;
        
        for i in range(len(self.comp)):
            endVal = endVal + self.comp[i].nParams
            if self.index[i] is not None:
                # only part of the data is involved in the kernel.
                g[startVal:endVal] = self.diagGradient(x[:, self.index[i]], covDiag=covGrad)
            else:
                # all the data is involved with the kernel.
                g[startVal:endVal] = self.diagGradient(x, covDiag=covGrad)
            startVal = endVal
        g = np.dot(g,self.paramGroups)        
        return g

    def display(self, numSpaces=0):
        
        """% DISPLAY Display parameters of the CMPND kernel.
        % FORMAT
        % DESC displays the parameters of the compound
        % kernel and the kernel type to the console.
        % ARG kern : the kernel to display.
        %
        % FORMAT does the same as above, but indents the display according
        % to the amount specified.
        % ARG kern : the kernel to display.
        % ARG spacing : how many spaces to indent the display of the kernel by.
        %
        % SEEALSO : cmpndKernParamInit, modelDisplay, kernDisplay
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """

        spacing = ' '*numSpaces
        print spacing, 'Compound kernel:'
        for i in range(len(self.comp)):
            self.comp[i].display(numSpaces+2)

    def gradX(self, X, X2):

        """% GRADX Gradient of CMPND kernel with respect to a point x.
        % FORMAT
        % DESC computes the gradient of the compound
        % kernel with respect to the input positions. 
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x : locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData x numInputs x numData. Where numData is
        % the number of data points and numInputs is the number of input
        % dimensions in X.
        %
        % FORMAT
        % DESC computes the gradident of the compound
        % kernel with respect to the input positions where both the row
        % positions and column positions are provided separately.
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x1 : row locations against which gradients are being computed.
        % ARG x2 : column locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData2 x numInputs x numData1. Where numData1 is
        % the number of data points in X1, numData2 is the number of data
        % points in X2 and numInputs is the number of input
        % dimensions in X.
        %
        % SEEALSO cmpndKernParamInit, kernGradX, cmpndKernDiagGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """
        

        i = 0

        if self.index[i] is not None:
            # only part of the data is involved with the kernel.
            gX = np.zeros((X2.shape[0], X2.shape[1], X.shape[0]))
            gX[:, self.index[i], :] = self.comp[i].gradX(X[:, self.index[i]], X2[:, self.index[i]])
        else:
            # all the data is involved with the kernel.
            gX = self.comp[i].gradX(X, X2)
        for i in range(1, len(self.comp)):
            if self.index[i] is not None:
                # only part of the data is involved with the kernel.
                gX[:, self.index[i], :] = gX[:, self.index[i], :] + self.comp[i].gradX(X[:, self.index[i]], X2[:, self.index[i]])
            else:
                # all the data is involved with the kernel.
                gX = gX + self.comp[i].gradX(X, X2)

        return gX

    def gradient(self, X, X2=None, covGrad=None):
        """% GRADIENT Gradient of CMPND kernel's parameters.
        % FORMAT
        % DESC computes the gradient of functions with respect to the
        % compound
        % kernel's parameters. As well as the kernel structure and the
        % input positions, the user provides a matrix PARTIAL which gives
        % the partial derivatives of the function with respect to the
        % relevant elements of the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x : the input locations for which the gradients are being
        % computed. 
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The argument takes
        % the form of a square matrix of dimension  numData, where numData is
        % the number of rows in X.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters. The ordering of the vector should match
        % that provided by the function kernExtractParam.
        %
        % FORMAT
        % DESC computes the derivatives as above, but input locations are
        % now provided in two matrices associated with rows and columns of
        % the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x1 : the input locations associated with the rows of the
        % kernel matrix.
        % ARG x2 : the input locations associated with the columns of the
        % kernel matrix.
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The matrix should
        % have the same number of rows as X1 and the same number of columns
        % as X2 has rows.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters.
        %
        % SEEALSO cmpndKernParamInit, kernGradient, cmpndKernDiagGradient, kernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """


        g = np.zeros(self.paramGroups.shape[0])
        startVal = 0
        endVal = 0

        for i in range(len(self.comp)):
            endVal = endVal + self.comp[i].nParams
            if self.index[i] is not None:
                # only part of the data is involved in the kernel.
                if X2 is None:
                    g[startVal:endVal]  = self.comp[i].gradient(X[:, self.index[i]], covGrad=covGrad)
                else:
                    g[startVal:endVal] = self.comp[i].gradient(X[:, self.index[i]], X2[:, self.index[i]], covGrad)
            else:
                # all the data is involved with the kernel.
                g[startVal:endVal]  = self.comp[i].gradient(X, X2, covGrad)
            startVal = endVal
        g = np.dot(g, self.paramGroups)
        return g
        

    def expandParam(self, params):
        """% EXPANDPARAM Create kernel structure from CMPND kernel's parameters.
        % FORMAT
        % DESC returns a compound kernel structure filled with the
        % parameters in the given vector. This is used as a helper function to
        % enable parameters to be optimised in, for example, the NETLAB
        % optimisation functions.
        % ARG kern : the kernel structure in which the parameters are to be
        % placed.
        % ARG param : vector of parameters which are to be placed in the
        % kernel structure.
        % RETURN kern : kernel structure with the given parameters in the
        % relevant locations.
        %
        % SEEALSO : cmpndKernParamInit, cmpndKernExtractParam, kernExpandParam
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """

        params = np.dot(params,self.paramGroups.T)
        startVal = 0
        endVal = 0
        self.whiteVariance = 0
        for i in range(len(self.comp)):
            endVal = endVal + self.comp[i].nParams
            self.comp[i].expandParam(params[startVal:endVal])
            startVal = endVal
            if self.comp[i].type[0:5] == 'white':
                # If kernel name starts with white, assume it is a white noise term.
                self.whiteVariance = self.whiteVariance + self.comp[i].variance
            else:
                if self.comp[i].whiteVariance is not None:
                    self.whiteVariance += self.comp[i].whiteVariance


    def extractParam(self):        
        """% EXTRACTPARAM Extract parameters from the CMPND kernel structure.
        % FORMAT
        % DESC Extract parameters from the compound kernel matrix into a vector of
        % parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the
        % kernel. The vector of 'transforms' is assumed to be empty
        % here. Any transformations of parameters should be done in
        % component kernels.
        %
        % SEEALSO cmpndKernParamInit, cmpndKernExpandParam, kernExtractParam, scg, conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """
        

        params = np.zeros(self.nParams);
        startVal = 0
        endVal = 0
        storedTypes = []
        for i in range(len(self.comp)):
            endVal = endVal + self.comp[i].nParams
            params[startVal:endVal] = self.comp[i].extractParam()
            startVal = endVal

        # If any parameters are 'tied together' deal with them.
        pGroups = self.paramGroups
        for i in range(self.paramGroups.shape[1]):
            ind = np.nonzero(pGroups[:, i])
            pGroups[ind[1:], i] = 0
        return np.dot(params,pGroups)

    def extractParamNames(self):
        
        """% EXTRACTPARAMNAMES Extract parameter names from the CMPND kernel structure.
        % FORMAT
        % DESC Extract parameters from the compound kernel matrix into a vector of
        % parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN names : a list of parameter names.
        %
        % SEEALSO cmpndKernParamInit, cmpndKernExpandParam, kernExtractParam, scg, conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """


        namesTemp = []
        names = []
        for i in range(self.nParams):
            namesTemp.append(None)
            names.append(None)

        startVal = 0
        endVal = 0
        storedTypes = []
        for i in range(len(self.comp)):
            endVal = endVal + self.comp[i].nParams
            namesTemp[startVal:endVal] = self.comp[i].extractParamNames()
            instNum = storedTypes.count(self.comp[i].type) + 1
            for j in range(startVal, endVal):
                namesTemp[j] = self.comp[i].type + ' ' + str(instNum) + ' ' + namesTemp[j]
            storedTypes.append(self.comp[i].type)
            startVal = endVal

        # If any parameters are 'tied together' deal with them.
        paramGroups = self.paramGroups
        for i in range(paramGroups.shape[1]):
            ind = np.nonzero(paramGroups[:, i])
            names[i] = namesTemp[ind[0]]
            if len(ind) > 1:
                for j in range(1,len(ind)):
                    names[i] = names[i] + ', ' + namesTemp[ind[j]]
            paramGroups[ind[1:], i] = 0
        return names

class rbf(kern):
    """% The radial basis function kernel (RBF) is sometimes also known as
    % the squared exponential kernel. It is a very smooth non-linear
    % kernel and is a popular choice for generic use.
    %
    % k(x_i, x_j) = sigma2 * exp(-gamma/2 *(x_i - x_j)'*(x_i - x_j))
    %
    % The parameters are sigma2, the process variance (kern.variance)
    % and gamma, the inverse width (kern.inverseWidth). The inverse
    % width controls how wide the basis functions are, the larger
    % gamma, the smaller the basis functions are.
    %
    % There is also an automatic relevance determination version of
    % this kernel provided.
    
    """
    
    def __init__(self, inDim=None, X=None):
        kern.__init__(self, inDim, X)

    def paramInit(self, inDim=None, X=None):
        """% RBFKERNPARAMINIT RBF kernel parameter initialisation.
        %
        % SEEALSO : rbfardKernParamInit
        %
        % FORMAT
        % DESC initialises the radial basis function
        %  kernel structure with some default parameters.
        % ARG kern : the kernel structure which requires initialisation.
        % RETURN kern : the kernel structure with the default parameters placed in.
        %
        % SEEALSO : create, kernParamInit
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
        
        """

        self.type = 'rbf'
        self.inverseWidth = 1.0
        self.variance = 1.0
        self.nParams = 2
        
        # Constrains parameters positive for optimisation.
        self.addTransform(optimi.defaultConstraint('positive'), [0, 1])
        self.stationary = True
        self.isNormalised = False

    def compute(self, x, x2=None):

        '''% DESC computes the kernel parameters for the radial basis function
        kernel given inputs associated with rows and columns.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : the input matrix associated with the rows of the kernel.
        % ARG x2 : the input matrix associated with the columns of the kernel.
        % RETURN k : the kernel matrix computed at the given points.
         
        % FORMAT
        % DESC computes the kernel matrix for the radial basis function
        kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % SEEALSO : diagCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009

        '''
        
        if x2 is None:
            n2 = dist2(x, x)
        else:
            n2 = dist2(x, x2)
        
        wi2 = (.5 * self.inverseWidth)
        sk = np.exp(-n2*wi2)
        k = sk*self.variance        
	if self.isNormalised:
	    k = k * np.sqrt(self.inverseWidth/(2.*np.pi))
        return k, sk, n2



    def diagCompute(self, x):
        '''% DIAGCOMPUTE Compute diagonal of RBF kernel.
        % FORMAT
        % DESC computes the diagonal of the kernel matrix for the radial basis function kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : a vector containing the diagonal of the kernel matrix
        % computed at the given points.
        %
        % SEEALSO : compute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
        
        '''

	if self.isNormalised:
		k = np.tile(self.variance * np.sqrt(self.inverseWidth/(2.*np.pi)), x.shape[0])
	else:
		k = np.tile(self.variance, x.shape[0])
	return k

    def diagGradX(self, X):
        '''% DIAGGRADX Gradient of RBF kernel's diagonal with respect to X.
        % FORMAT
        % DESC computes the gradient of the diagonal of the radial basis function kernel matrix with
        % respect to the elements of the design matrix given in X.
        % ARG kern : the kernel structure for which gradients are being computed.
        % ARG X : the input data in the form of a design matrix.
        % RETURN gX : the gradients of the diagonal with respect to each element
        % of X. The returned matrix has the same dimensions as X.
        %
        % SEEALSO : gradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009

        '''
        
        return np.zeros(X.shape)


    def diagGradient(self, X, covDiag):
        '''% DIAGGRADIENT Compute the gradient of the RBF kernel's diagonal wrt parameters.
        % FORMAT
        % DESC computes the gradient of functions of the diagonal of the
        % radial basis function kernel matrix with respect to the parameters of the kernel. The
        % parameters' gradients are returned in the order given by the
        % rbfKernExtractParam command.
        % ARG kern : the kernel structure for which the gradients are
        % computed.
        % ARG x : the input data for which the gradient is being computed.
        % ARG factors : partial derivatives of the function of interest with
        % respect to the diagonal elements of the kernel.
        % RETURN g : gradients of the relevant function with respect to each
        % of the parameters. Ordering should match the ordering given in
        % rbfKernExtractParam.
        %
        % SEEALSO : gradient
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009

        '''

        g = np.zeros(self.nParams)
	if self.isNormalised:
		g[0] = np.sqrt(self.inverseWidth/(2*np.pi)) * covDiag.sum()
		g[1] = 0.5 * self.variance / np.sqrt(2*np.pi*self.inverseWidth) * covDiag.sum()
	else:
		g[0] = 0.0
		g[1] = covDiag.sum()
        return g

    def display(self, numSpaces=0):
        '''% DISPLAY Display parameters of the RBF kernel.
        % FORMAT
        % DESC displays the parameters of the radial basis function
        % kernel and the kernel type to the console.
        % ARG kern : the kernel to display.
        %
        % FORMAT does the same as above, but indents the display according
        % to the amount specified.
        % ARG kern : the kernel to display.
        % ARG spacing : how many spaces to indent the display of the kernel by.
        %
        % SEEALSO :
        %
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009

        '''
        
        spacing = ' '*numSpaces
        print spacing, 'RBF inverse width: ', self.inverseWidth, '(length scale ', 1.0/math.sqrt(self.inverseWidth), ')'
        print spacing, 'RBF variance: ', self.variance

    def expandParam(self, params):
        '''% EXPANDPARAM Create kernel structure from RBF kernel's parameters.
        % FORMAT
        % DESC returns a radial basis function kernel structure filled with the
        % parameters in the given vector. This is used as a helper function to
        % enable parameters to be optimised in, for example, the NETLAB
        % optimisation functions.
        % ARG kern : the kernel structure in which the parameters are to be
        % placed.
        % ARG param : vector of parameters which are to be placed in the
        % kernel structure.
        % RETURN kern : kernel structure with the given parameters in the
        % relevant locations.
        %
        % SEEALSO : extractParam
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009

        '''

        self.inverseWidth = params[0]
        self.variance = params[1]

    def extractParam(self):
        '''% EXTRACTPARAM Extract parameters from the RBF kernel structure.
        % FORMAT
        % DESC Extract parameters from the radial basis function kernel
        % structure into a vector of parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        %
        % SEEALSO expandParam, netlab.scg, netlab.conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009

        '''
        
        return np.array([self.inverseWidth, self.variance])

    def extractParamNames(self):
        '''% EXTRACTPARAMNAMES Extract parameter names from the RBF kernel structure.
        % FORMAT
        % DESC Extract parameter names from the radial basis
        % function kernel structure.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN names : cell array of strings giving names to the parameters.
        
        '''
        
        return ['inverse width', 'variance']

    def gradX(self, X, X2):
        '''% GRADX Gradient of RBF kernel with respect to input locations.
        % FORMAT
        % DESC computes the gradident of the radial basis function
        % kernel with respect to the input positions where both the row
        % positions and column positions are provided separately.
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x1 : row locations against which gradients are being computed.
        % ARG x2 : column locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData2 x numInputs x numData1. Where numData1 is
        % the number of data points in X1, numData2 is the number of data
        % points in X2 and numInputs is the number of input
        % dimensions in X.
        %
        % SEEALSO : diagGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009

        '''
        
        gX = np.empty((X2.shape[0], X2.shape[1], X.shape[0]))
        for i in range(X.shape[0]):
            gX[:, :, i] = self.gradXpoint(X[i:i+1, :], X2)
	if self.isNormalised:
		gX = gX * np.sqrt(self.inverseWidth/(2*np.pi))
        return gX


    def gradXpoint(self, x, X2):
        '''% GRADXPOINT Gradient with respect to one point of x.'''

        gX = np.zeros(X2.shape)
        n2 = dist2(X2, x)
        wi2 = (.5 * self.inverseWidth)
        rbfPart = self.variance*np.exp(-n2*wi2)
        for i in range(x.shape[1]):
            gX[:, i:i+1] = self.inverseWidth*(X2[:, i:i+1] - x[:,i])*rbfPart
        return gX

    def gradient(self, X, X2=None, covGrad=None):
        '''% GRADIENT Gradient of RBF kernel's parameters.
        % FORMAT
        % DESC computes the gradient of functions with respect to the
        % radial basis function
        % kernel's parameters. As well as the kernel structure and the
        % input positions, the user provides a matrix PARTIAL which gives
        % the partial derivatives of the function with respect to the
        % relevant elements of the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x : the input locations for which the gradients are being
        % computed. 
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The argument takes
        % the form of a square matrix of dimension  numData, where numData is
        % the number of rows in X.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters. The ordering of the vector should match
        % that provided by the function extractParam.
        %
        % FORMAT
        % DESC computes the derivatives as above, but input locations are
        % now provided in two matrices associated with rows and columns of
        % the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x1 : the input locations associated with the rows of the
        % kernel matrix.
        % ARG x2 : the input locations associated with the columns of the
        % kernel matrix.
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The matrix should
        % have the same number of rows as X1 and the same number of columns
        % as X2 has rows.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters.
        %
        % SEEALSO diagGradient, gradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009

        '''
        if X2 is None:
            k, sk, dist2xx = self.compute(X)
        else:
            k, sk, dist2xx = self.compute(X, X2)
        g = np.zeros(self.nParams)
	g[0] = - .5*np.multiply(np.multiply(covGrad,k),dist2xx).sum()
	g[1] =  np.multiply(covGrad,sk).sum()
	if self.isNormalised:
		g[0] = (0.5*self.variance/np.sqrt(2*np.pi)) * np.sum(covgrad*sk* (1./np.sqrt(self.inverseWidth)-np.sqrt(self.inverseWidth)*dist2xx))
		g[1] *= np.sqrt(self.inverseWidth/(2*np.pi)) 
        return g

class lin(kern):
    """% The linear kernel (LIN) is the simple inner product
    % kernel. Sampling from this kernel produces linear functions.
    %
    % k(x_i, x_j) = sigma2 * x_i'*x_j
    %
    % There is one parameter, sigma2, which is stored in the field
    % kern.variance.
    
    """
    
    def __init__(self, inDim=None, X=None):
        kern.__init__(self, inDim, X)

    def paramInit(self, inDim=None, X=None):
        """% LINKERNPARAMINIT LIN kernel parameter initialisation.
        % The linear kernel (LIN) is the simple inner product
        % kernel. Sampling from this kernel produces linear functions.
        %
        % k(x_i, x_j) = sigma2 * x_i'*x_j
        %
        % There is one parameter, sigma2, which is stored in the field
        % kern.variance.
        %
        % SEEALSO : linardKernParamInit
        %
        % FORMAT
        % DESC initialises the linear
        %  kernel structure with some default parameters.
        % ARG kern : the kernel structure which requires initialisation.
        % RETURN kern : the kernel structure with the default parameters placed in.
        %
        % SEEALSO : kernCreate, kernParamInit
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """

        self.type = 'lin'
        self.variance = 1.0
        self.nParams = 1
        
        # Constrains parameters positive for optimisation.
        self.addTransform(optimi.defaultConstraint('positive'), [0])
        self.stationary = False
        self.normalised = False

    def compute(self, x, x2=None):

        """% LINKERNCOMPUTE Compute the LIN kernel given the parameters and X.
        % FORMAT
        % DESC computes the kernel parameters for the linear
        % kernel given inputs associated with rows and columns.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : the input matrix associated with the rows of the kernel.
        % ARG x2 : the input matrix associated with the columns of the kernel.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % FORMAT
        % DESC computes the kernel matrix for the linear
        % kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % SEEALSO : linKernParamInit, kernCompute, kernCreate, linKernDiagCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """


        if x2 is None:
            sk = np.dot(x,x.T) 
        else:
            sk = np.dot(x,x2.T)
        k = sk*self.variance
        return k, sk



    def diagCompute(self, x):
        """% LINKERNDIAGCOMPUTE Compute diagonal of LIN kernel.
        % FORMAT
        % DESC computes the diagonal of the kernel matrix for the linear kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : a vector containing the diagonal of the kernel matrix
        % computed at the given points.
        %
        % SEEALSO : linKernParamInit, kernDiagCompute, kernCreate, linKernCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return np.sum(x*x, 1)*self.variance

    def diagGradX(self, X):

        """% LINKERNDIAGGRADX Gradient of LIN kernel's diagonal with respect to X.
        % FORMAT
        % DESC computes the gradient of the diagonal of the linear kernel matrix with
        % respect to the elements of the design matrix given in X.
        % ARG kern : the kernel structure for which gradients are being computed.
        % ARG X : the input data in the form of a design matrix.
        % RETURN gX : the gradients of the diagonal with respect to each element
        % of X. The returned matrix has the same dimensions as X.
        %
        % SEEALSO : linKernParamInit, kernDiagGradX, linkernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return 2*X*self.variance


    def diagGradient(self, X, covDiag):
        '''% DIAGGRADIENT Compute the gradient of the RBF kernel's diagonal wrt parameters.
        % FORMAT
        % DESC computes the gradient of functions of the diagonal of the
        % radial basis function kernel matrix with respect to the parameters of the kernel. The
        % parameters' gradients are returned in the order given by the
        % rbfKernExtractParam command.
        % ARG kern : the kernel structure for which the gradients are
        % computed.
        % ARG x : the input data for which the gradient is being computed.
        % ARG factors : partial derivatives of the function of interest with
        % respect to the diagonal elements of the kernel.
        % RETURN g : gradients of the relevant function with respect to each
        % of the parameters. Ordering should match the ordering given in
        % rbfKernExtractParam.
        %
        % SEEALSO : gradient
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''


        g = np.zeros(self.nParams)
        g[0] = 0.0
        g[1] = covDiag.sum()
        return g

    def display(self, numSpaces=0):
        '''% DISPLAY Display parameters of the LIN kernel.
        % FORMAT
        % DESC displays the parameters of the radial basis function
        % kernel and the kernel type to the console.
        % ARG kern : the kernel to display.
        %
        % FORMAT does the same as above, but indents the display according
        % to the amount specified.
        % ARG kern : the kernel to display.
        % ARG spacing : how many spaces to indent the display of the kernel by.
        %
        % SEEALSO :
        %
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009
        
        '''
        
        spacing = ' '*numSpaces
        print spacing, 'Linear variance: ', self.variance

    def expandParam(self, params):
        '''% EXPANDPARAM Create kernel structure from LIN kernel's parameters.
        % FORMAT
        % DESC returns a linear kernel structure filled with the
        % parameters in the given vector. This is used as a helper function to
        % enable parameters to be optimised in, for example, the NETLAB
        % optimisation functions.
        % ARG kern : the kernel structure in which the parameters are to be
        % placed.
        % ARG param : vector of parameters which are to be placed in the
        % kernel structure.
        % RETURN kern : kernel structure with the given parameters in the
        % relevant locations.
        %
        % SEEALSO : extractParam
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''

        self.variance = params[0]

    def extractParam(self):
        '''% EXTRACTPARAM Extract parameters from the LIN kernel structure.
        % FORMAT
        % DESC Extract parameters from the linear kernel
        % structure into a vector of parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        %
        % SEEALSO expandParam, netlab.scg, netlab.conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009

        '''
        
        return np.array([self.variance])

    def extractParamNames(self):
        '''% EXTRACTPARAMNAMES Extract parameter names from the LIN kernel structure.
        % FORMAT
        % DESC Extract parameter names from the linear kernel structure.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN names : cell array of strings giving names to the parameters.'''
        return ['variance']

    def gradX(self, X, X2):
        '''% GRADX Gradient of LIN kernel with respect to input locations.
        % FORMAT
        % DESC computes the gradient of the linear
        % kernel with respect to the input positions where both the row
        % positions and column positions are provided separately.
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x1 : row locations against which gradients are being computed.
        % ARG x2 : column locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData2 x numInputs x numData1. Where numData1 is
        % the number of data points in X1, numData2 is the number of data
        % points in X2 and numInputs is the number of input
        % dimensions in X.
        %
        % SEEALSO : diagGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''
        
        
        return np.tile(np.reshape(self.variance*X2, (X2.shape[0], X2.shape[1], 1)), (1, 1, X.shape[0]))



    def gradient(self, X, X2=None, covGrad=None):

        """% LINKERNGRADIENT Gradient of LIN kernel's parameters.
        % FORMAT
        % DESC computes the gradient of functions with respect to the
        % linear
        % kernel's parameters. As well as the kernel structure and the
        % input positions, the user provides a matrix PARTIAL which gives
        % the partial derivatives of the function with respect to the
        % relevant elements of the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x : the input locations for which the gradients are being
        % computed. 
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The argument takes
        % the form of a square matrix of dimension  numData, where numData is
        % the number of rows in X.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters. The ordering of the vector should match
        % that provided by the function kernExtractParam.
        %
        % FORMAT
        % DESC computes the derivatives as above, but input locations are
        % now provided in two matrices associated with rows and columns of
        % the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x1 : the input locations associated with the rows of the
        % kernel matrix.
        % ARG x2 : the input locations associated with the columns of the
        % kernel matrix.
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The matrix should
        % have the same number of rows as X1 and the same number of columns
        % as X2 has rows.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters.
        %
        % SEEALSO linKernParamInit, kernGradient, linKernDiagGradient, kernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """

        if X2 is None:
            k, sk = self.compute(X)
        else:
            k, sk = self.compute(X, X2)
        g = np.zeros(self.nParams)
        g[0] =  np.multiply(covGrad,sk).sum()
        return g

    
class white(kern):
    """% The white noise kernel arises from assuming independent Gaussian
    % noise for each point in the function. The variance of the noise is
    % given by the kern.variance parameter.
    % 
    % This kernel is not intended to be used independently, it is provided
    % so that it may be combined with other kernels in a compound kernel."""

    def __init__(self, inDim=None, X=None):
        kern.__init__(self, inDim, X)

    def paramInit(self, inDim=None, X=None):
        """% WHITEKERNPARAMINIT WHITE kernel parameter initialisation.
        %
        % SEEALSO : cmpndKernParamInit
        %
        % FORMAT
        % DESC initialises the white noise
        %  kernel structure with some default parameters.
        % ARG kern : the kernel structure which requires initialisation.
        % RETURN kern : the kernel structure with the default parameters placed in.
        %
        % SEEALSO : create, kernParamInit
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        self.type = 'white'
        self.variance = math.exp(-2.0)
        self.nParams = 1
        
        # Constrains parameters positive for optimisation.
        self.addTransform(optimi.defaultConstraint('positive'), [0])
        self.stationary = True

    def compute(self, x, x2=None):
        """% WHITEKERNCOMPUTE Compute the WHITE kernel given the parameters and X.
        % FORMAT
        % DESC computes the kernel parameters for the white noise
        % kernel given inputs associated with rows and columns.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : the input matrix associated with the rows of the kernel.
        % ARG x2 : the inpute matrix associated with the columns of the kernel.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % FORMAT
        % DESC computes the kernel matrix for the white noise
        % kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % SEEALSO : whiteKernParamInit, kernCompute, create, whiteKernDiagCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        if x2 is None:
            sk = np.eye(x.shape[0])
            k = self.variance*sk;
        else:
            sk = np.zeros((x.shape[0], x2.shape[0]))
            k = np.zeros((x.shape[0], x2.shape[0]))

        return k, sk

    def diagCompute(self, x):  
        """% WHITEKERNDIAGCOMPUTE Compute diagonal of WHITE kernel.
        % FORMAT
        % DESC computes the diagonal of the kernel matrix for the white noise kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : a vector containing the diagonal of the kernel matrix
        % computed at the given points.
        %
        % SEEALSO : whiteKernParamInit, kernDiagCompute, create, whiteKernCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """
        
        return np.tile(self.variance, x.shape[0])


    def diagGradX(self, X):
        """% WHITEKERNDIAGGRADX Gradient of WHITE kernel's diagonal with respect to X.
        % FORMAT
        % DESC computes the gradient of the diagonal of the white noise kernel matrix with
        % respect to the elements of the design matrix given in X.
        % ARG kern : the kernel structure for which gradients are being computed.
        % ARG X : the input data in the form of a design matrix.
        % RETURN gX : the gradients of the diagonal with respect to each element
        % of X. The returned matrix has the same dimensions as X.
        %
        % SEEALSO : whiteKernParamInit, kernDiagGradX, whitekernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return np.zeros(X.shape)

    def diagGradient(self, X, covDiag):
        """% WHITEKERNDIAGGRADIENT Compute the gradient of the WHITE kernel's diagonal wrt parameters.
        % FORMAT
        % DESC computes the gradient of functions of the diagonal of the
        % white noise kernel matrix with respect to the parameters of the kernel. The
        % parameters' gradients are returned in the order given by the
        % whiteKernExtractParam command.
        % ARG kern : the kernel structure for which the gradients are
        % computed.
        % ARG x : the input data for which the gradient is being computed.
        % ARG factors : partial derivatives of the function of interest with
        % respect to the diagonal elements of the kernel.
        % RETURN g : gradients of the relevant function with respect to each
        % of the parameters. Ordering should match the ordering given in
        % whiteKernExtractParam.
        %
        % SEEALSO : whiteKernParamInit, kernDiagGradient, whiteKernExtractParam, whiteKernGradient
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return np.array([np.sum(covDiag)])

    def display(self, numSpaces=0):
        """% WHITEKERNDISPLAY Display parameters of the WHITE kernel.
        % FORMAT
        % DESC displays the parameters of the white noise
        % kernel and the kernel type to the console.
        % ARG kern : the kernel to display.
        %
        % FORMAT does the same as above, but indents the display according
        % to the amount specified.
        % ARG kern : the kernel to display.
        % ARG spacing : how many spaces to indent the display of the kernel by.
        %
        % SEEALSO : whiteKernParamInit, modelDisplay, kernDisplay
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        spacing = ' '*numSpaces
        print spacing, 'White Noise Variance: ', self.variance

    def expandParam(self, params):
        """% WHITEKERNEXPANDPARAM Create kernel structure from WHITE kernel's parameters.
        % FORMAT
        % DESC returns a white noise kernel structure filled with the
        % parameters in the given vector. This is used as a helper function to
        % enable parameters to be optimised in, for example, the NETLAB
        % optimisation functions.
        % ARG kern : the kernel structure in which the parameters are to be
        % placed.
        % ARG param : vector of parameters which are to be placed in the
        % kernel structure.
        % RETURN kern : kernel structure with the given parameters in the
        % relevant locations.
        %
        % SEEALSO : whiteKernParamInit, whiteKernExtractParam, kernExpandParam
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
        
        """

        self.variance = params[0]

    def extractParam(self):
        """% WHITEKERNEXTRACTPARAM Extract parameters from the WHITE kernel structure.
        % FORMAT
        % DESC Extract parameters from the white noise kernel structure into a
        % vector of parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        %
        % FORMAT
        % DESC Extract parameters and parameter names from the white noise
        % kernel structure.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        % RETURN names : cell array of strings giving paramter names.
        %
        % SEEALSO whiteKernParamInit, whiteKernExpandParam, kernExtractParam, scg, conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
        %
        """
        
        return np.array([self.variance])

    def extractParamNames(self):

        return ['variance']

    def gradX(self, X, X2=None):
        
        """% WHITEKERNGRADX Gradient of WHITE kernel with respect to input locations.
        % FORMAT
        % DESC computes the gradident of the white noise
        % kernel with respect to the input positions where both the row
        % positions and column positions are provided separately.
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x1 : row locations against which gradients are being computed.
        % ARG x2 : column locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData2 x numInputs x numData1. Where numData1 is
        % the number of data points in X1, numData2 is the number of data
        % points in X2 and numInputs is the number of input
        % dimensions in X.
        %
        % SEEALSO whiteKernParamInit, kernGradX, whiteKernDiagGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
        
        """

        if X2 is None:
            return np.zeros((X.shape[0], X.shape[1], X.shape[0]))
        else:
            return np.zeros((X2.shape[0], X2.shape[1], X.shape[0]))

    def gradient(self, X, X2=None, covGrad=None):

        """% WHITEKERNGRADIENT Gradient of WHITE kernel's parameters.
        % FORMAT
        % DESC computes the gradient of functions with respect to the
        % white noise
        % kernel's parameters. As well as the kernel structure and the
        % input positions, the user provides a matrix PARTIAL which gives
        % the partial derivatives of the function with respect to the
        % relevant elements of the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x : the input locations for which the gradients are being
        % computed. 
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The argument takes
        % the form of a square matrix of dimension  numData, where numData is
        % the number of rows in X.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters. The ordering of the vector should match
        % that provided by the function kernExtractParam.
        %
        % FORMAT
        % DESC computes the derivatives as above, but input locations are
        % now provided in two matrices associated with rows and columns of
        % the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x1 : the input locations associated with the rows of the
        % kernel matrix.
        % ARG x2 : the input locations associated with the columns of the
        % kernel matrix.
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The matrix should
        % have the same number of rows as X1 and the same number of columns
        % as X2 has rows.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters.
        %
        % SEEALSO whiteKernParamInit, kernGradient, whiteKernDiagGradient, kernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """


        if X2 is None:
            return np.array([np.trace(covGrad)])
        else:
            return np.array([0.0])



class bias(kern):
    """% The white noise kernel arises from assuming independent Gaussian
    % noise for each point in the function. The variance of the noise is
    % given by the kern.variance parameter.
    % 
    % This kernel is not intended to be used independently, it is provided
    % so that it may be combined with other kernels in a compound kernel."""

    def __init__(self, inDim=None, X=None):
        kern.__init__(self, inDim, X)

    def paramInit(self, inDim=None, X=None):
        """% BIASKERNPARAMINIT BIAS kernel parameter initialisation.
        %
        % SEEALSO : cmpndKernParamInit
        %
        % FORMAT
        % DESC initialises the bias noise
        %  kernel structure with some default parameters.
        % ARG kern : the kernel structure which requires initialisation.
        % RETURN kern : the kernel structure with the default parameters placed in.
        %
        % SEEALSO : create, kernParamInit
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """
        self.type = 'bias'
        self.variance = math.exp(-2.0)
        self.nParams = 1
        
        # Constrains parameters positive for optimisation.
        self.addTrans(optimi.defaultConstraint('positive'), [0])
        self.stationary = True

    def compute(self, x, x2=None):
        """% BIASKERNCOMPUTE Compute the BIAS kernel given the parameters and X.
        % FORMAT
        % DESC computes the kernel parameters for the bias noise
        % kernel given inputs associated with rows and columns.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : the input matrix associated with the rows of the kernel.
        % ARG x2 : the inpute matrix associated with the columns of the kernel.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % FORMAT
        % DESC computes the kernel matrix for the bias noise
        % kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : the kernel matrix computed at the given points.
        %
        % SEEALSO : biasKernParamInit, kernCompute, create, biasKernDiagCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        if x2 is None:
            sk = np.ones((x.shape[0], x.shape[0]))
            k = np.tile(self.variance, (x.shape[0], x.shape[0]))
        else:
            sk = np.ones((x.shape[0], x2.shape[0]))
            k = np.tile(self.variance, (x.shape[0], x2.shape[0]))

        return k, sk

    def diagCompute(self, x):  
        """% BIASKERNDIAGCOMPUTE Compute diagonal of BIAS kernel.
        % FORMAT
        % DESC computes the diagonal of the kernel matrix for the bias noise kernel given a design matrix of inputs.
        % ARG kern : the kernel structure for which the matrix is computed.
        % ARG x : input data matrix in the form of a design matrix.
        % RETURN k : a vector containing the diagonal of the kernel matrix
        % computed at the given points.
        %
        % SEEALSO : biasKernParamInit, kernDiagCompute, create, biasKernCompute
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """
        
        return np.tile(self.variance, x.shape[0])


    def diagGradX(self, X):
        """% BIASKERNDIAGGRADX Gradient of BIAS kernel's diagonal with respect to X.
        % FORMAT
        % DESC computes the gradient of the diagonal of the bias noise kernel matrix with
        % respect to the elements of the design matrix given in X.
        % ARG kern : the kernel structure for which gradients are being computed.
        % ARG X : the input data in the form of a design matrix.
        % RETURN gX : the gradients of the diagonal with respect to each element
        % of X. The returned matrix has the same dimensions as X.
        %
        % SEEALSO : biasKernParamInit, kernDiagGradX, biaskernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return np.zeros(X.shape)

    def diagGradient(self, X, covDiag):
        """% BIASKERNDIAGGRADIENT Compute the gradient of the BIAS kernel's diagonal wrt parameters.
        % FORMAT
        % DESC computes the gradient of functions of the diagonal of the
        % bias noise kernel matrix with respect to the parameters of the kernel. The
        % parameters' gradients are returned in the order given by the
        % biasKernExtractParam command.
        % ARG kern : the kernel structure for which the gradients are
        % computed.
        % ARG x : the input data for which the gradient is being computed.
        % ARG factors : partial derivatives of the function of interest with
        % respect to the diagonal elements of the kernel.
        % RETURN g : gradients of the relevant function with respect to each
        % of the parameters. Ordering should match the ordering given in
        % biasKernExtractParam.
        %
        % SEEALSO : biasKernParamInit, kernDiagGradient, biasKernExtractParam, biasKernGradient
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        return np.array([np.sum(covDiag)])

    def display(self, numSpaces=0):
        """% BIASKERNDISPLAY Display parameters of the BIAS kernel.
        % FORMAT
        % DESC displays the parameters of the bias noise
        % kernel and the kernel type to the console.
        % ARG kern : the kernel to display.
        %
        % FORMAT does the same as above, but indents the display according
        % to the amount specified.
        % ARG kern : the kernel to display.
        % ARG spacing : how many spaces to indent the display of the kernel by.
        %
        % SEEALSO : biasKernParamInit, modelDisplay, kernDisplay
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
        
        """

        spacing = ' '*numSpaces
        print spacing, 'Bias Variance: ', self.variance

    def expandParam(self, params):
        """% BIASKERNEXPANDPARAM Create kernel structure from BIAS kernel's parameters.
        % FORMAT
        % DESC returns a bias noise kernel structure filled with the
        % parameters in the given vector. This is used as a helper function to
        % enable parameters to be optimised in, for example, the NETLAB
        % optimisation functions.
        % ARG kern : the kernel structure in which the parameters are to be
        % placed.
        % ARG param : vector of parameters which are to be placed in the
        % kernel structure.
        % RETURN kern : kernel structure with the given parameters in the
        % relevant locations.
        %
        % SEEALSO : biasKernParamInit, biasKernExtractParam, kernExpandParam
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
                
        """

        self.variance = params[0]

    def extractParam(self):
        """% BIASKERNEXTRACTPARAM Extract parameters from the BIAS kernel structure.
        % FORMAT
        % DESC Extract parameters from the bias noise kernel structure into a
        % vector of parameters for optimisation.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        %
        % FORMAT
        % DESC Extract parameters and parameter names from the bias noise
        % kernel structure.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN param : vector of parameters extracted from the kernel. If
        % the field 'transforms' is not empty in the kernel matrix, the
        % parameters will be transformed before optimisation (for example
        % positive only parameters could be logged before being returned).
        % RETURN names : cell array of strings giving paramter names.
        %
        % SEEALSO biasKernParamInit, biasKernExpandParam, kernExtractParam, scg, conjgrad
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

        """
        
        return np.array([self.variance])

    def extractParamNames(self):

        return ['variance']

    def gradX(self, X, X2=None):
        
        """% BIASKERNGRADX Gradient of BIAS kernel with respect to input locations.
        % FORMAT
        % DESC computes the gradident of the bias noise
        % kernel with respect to the input positions where both the row
        % positions and column positions are provided separately.
        % ARG kern : kernel structure for which gradients are being
        % computed.
        % ARG x1 : row locations against which gradients are being computed.
        % ARG x2 : column locations against which gradients are being computed.
        % RETURN g : the returned gradients. The gradients are returned in
        % a matrix which is numData2 x numInputs x numData1. Where numData1 is
        % the number of data points in X1, numData2 is the number of data
        % points in X2 and numInputs is the number of input
        % dimensions in X.
        %
        % SEEALSO biasKernParamInit, kernGradX, biasKernDiagGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

        """

        if X2 is None:
            return np.zeros((X.shape[0], X.shape[1], X.shape[0]))
        else:
            return np.zeros((X2.shape[0], X2.shape[1], X.shape[0]))

    def gradient(self, X, X2=None, covGrad=None):

        """% BIASKERNGRADIENT Gradient of BIAS kernel's parameters.
        % FORMAT
        % DESC computes the gradient of functions with respect to the
        % bias noise
        % kernel's parameters. As well as the kernel structure and the
        % input positions, the user provides a matrix PARTIAL which gives
        % the partial derivatives of the function with respect to the
        % relevant elements of the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x : the input locations for which the gradients are being
        % computed. 
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The argument takes
        % the form of a square matrix of dimension  numData, where numData is
        % the number of rows in X.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters. The ordering of the vector should match
        % that provided by the function kernExtractParam.
        %
        % FORMAT
        % DESC computes the derivatives as above, but input locations are
        % now provided in two matrices associated with rows and columns of
        % the kernel matrix. 
        % ARG kern : the kernel structure for which the gradients are being
        % computed.
        % ARG x1 : the input locations associated with the rows of the
        % kernel matrix.
        % ARG x2 : the input locations associated with the columns of the
        % kernel matrix.
        % ARG partial : matrix of partial derivatives of the function of
        % interest with respect to the kernel matrix. The matrix should
        % have the same number of rows as X1 and the same number of columns
        % as X2 has rows.
        % RETURN g : gradients of the function of interest with respect to
        % the kernel parameters.
        %
        % SEEALSO biasKernParamInit, kernGradient, biasKernDiagGradient, kernGradX
        %
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

        """

        return np.array([np.sum(covGrad)])

class sim(kern):
	def __init__(self, inDim=None, X=None):
		kern.__init__(self,inDim, X)

	def paramInit(self,inDim):
		"""% SIMKERNPARAMINIT SIM kernel parameter initialisation.
		% The single input motif (SIM) kernel is specifically designed for
		% working with gene networks where there is assumed to be a single
		% transcription factor controlling several genes. If each gene is
		% related to the transcription factor through the following
		% differential equation,
		%
		% dx(t)/dt = B + S f(t-delta) - D x(t),
		%
		% where D is a decay term, S is a response term, delta is a time delay
		% and B is an initial level. Then if f(t) is assumed to come from a
		% Gaussian process with an RBF covariance function x(t) is a Gaussian
		% process with a covariance function provided by the single input
		% motif kernel.
		%
		% The kernel is designed to interoperate with the multiple output
		% block kernel so that f(t) can be inferred given several different
		% instantiations of x(t) (associated with different genes).
		%
		% By default the parameters (B, S, delta and D) are constrained positive. If
		% kern.options.isNegativeS is set true then the parameter S is allowed to go
		% negative.
		%
		% SEEALSO : kernCreate, kernParamInit, simKernCompute
		%
		% COPYRIGHT : Neil D. Lawrence, 2006, 2009

		"""
		assert(self.inputDimension==1,'SIM kernel only valid for one-D input.')

		self.gaussianInitial = False;
		if hasattr(self, 'options'):
			if hasattr(self.options, 'gaussianInitial'):
				if self.options.gaussianInitial:
					self.gaussianInitial = True;
					self.initialVariance = 1;

		self.delay = 0.
		self.decay = 1.
		self.initVal = 1.
		self.variance = 1.
		self.inverseWidth = 1.

		if self.gaussianInitial:
			self.nParams = 4
		else:
			self.nParams = 3

		self.isNegativeS = False
		if hasattr(self, 'options'):
			if hasattr(self.options, 'isNegativeS'):
				if self.options.isNegativeS:
					self.isNegativeS = True

		if self.isNegativeS:
			index = np.setdiff1d(np.arange(self.nParams), np.array([2]))
			self.addTransform(optimi.defaultConstraint('positive'), index)
			self.sensitivity = 1.
		else:
			self.addTransform(optimi.defaultConstraint('positive'), np.arange(self.nParams))

		self.isStationary = False;
		self.isNormalised = False;
		self.positiveTime = True;

	def compute(self,t,t2=None):
		"""% SIMKERNCOMPUTE Compute the SIM kernel given the parameters and X.
		% FORMAT
		% DESC computes the kernel parameters for the single input motif
		% kernel given inputs associated with rows and columns.
		% ARG kern : the kernel structure for which the matrix is computed.
		% ARG t1 : the input matrix associated with the rows of the kernel.
		% ARG t2 : the input matrix associated with the columns of the kernel.
		% RETURN k : the kernel matrix computed at the given points.
		% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
		% part).
		%
		% FORMAT
		% DESC computes the kernel matrix for the single input motif
		% kernel given a design matrix of inputs.
		% ARG kern : the kernel structure for which the matrix is computed.
		% ARG t : input data matrix in the form of a design matrix.
		% RETURN k : the kernel matrix computed at the given points.
		% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
		% part).
		%
		% SEEALSO : simKernParamInit, kernCompute, kernCreate, simKernDiagCompute
		%
		% COPYRIGHT : Neil D. Lawrence, 2006
		%
		% MODIFICATIONS : David Luengo, 2009
		%
		% MODIFICATIONS : Mauricio Alvarez, 2009
		%
		% MODIFICATIONS : Antti Honkela, 2009
		%
		% MODIFICATIONS : James Hensman 2011
		"""
		symmetry=False
		if t2 is None:
			t2 = t
			symmetry=True


		assert( (t.shape[1]==1) & (t2.shape[1]==1),'Input can only have one column')

		sigma = np.sqrt(2/self.inverseWidth)

		if (self.isStationary):
			h = simComputeHStat(t, t2, self.decay, self.decay, self.delay, self.delay, sigma)[0]
		else:
			h = simComputeH(t, t2, self.decay, self.decay, self.delay, self.delay, sigma)[0]
		if symmetry:
			sk = 0.5 * (h + h.T)
		else:
			if (self.isStationary == False):
				h2,tmp1,tmp2,tmp3 = simComputeH(t2, t, self.decay, self.decay, self.delay, self.delay, sigma)
			else:
				h2,tmp1,tmp2,tmp3 = simComputeHStat(t2, t, self.decay, self.decay, self.delay, self.delay, sigma)
			sk = 0.5 * (h + h2.T)

		#if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
		if self.isNormalised:
			k = sk
		else:
			k = np.sqrt(np.pi)*sigma*sk

		#if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
		if self.isNegativeS:
			k = (kern.sensitivity*kern.sensitivity)*k
		else:
			k = self.variance*k

		#if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
		if self.gaussianInitial:
			dim1 = t.shape[0]
			dim2 = t2.shape[0]
			t1Mat = t*np.ones(1,dim2)#(:, ones(1, dim2))
			t2Mat = t2*np.ones(1,dim1)#(:, ones(1, dim1))'
			t2Mat = t2Mat.T

			k = k + self.initialVariance * exp(- self.decay * (t1Mat + t2Mat))
		return(k)

	def expandParam(self,param):
		self.decay,self.inverseWidth, self.variance = param

	def extractParam(self):
		return np.array([self.decay,self.inverseWidth, self.variance])

	def crossCompute(self,other,t1,t2=None):
		if (t2 is None):
			t2=t1

		if isinstance(other, sim):
			assert self.isNormalised == other.isNormalised
			assert self.inverseWidth == other.inverseWidth

			sigma = np.sqrt(2./self.inverseWidth)

			if not self.isStationary:
				h1,tmp1,tmp2,tmp3 = simComputeH(t1, t2, self.decay, other.decay, self.delay, other.delay, sigma)
			        h2,tmp1,tmp2,tmp3 = simComputeH(t2, t1, other.decay, self.decay, other.delay, self.delay, sigma);
			else:
				h1,tmp1,tmp2,tmp3 = simComputeHStat(t1, t2, self.decay, other.decay, self.delay, other.delay, sigma);
			        h2,tmp1,tmp2,tmp3 = simComputeHStat(t2, t1, other.decay, self.decay, other.delay, self.delay, sigma);

			sK = 0.5 * (h1 + h2.T)

			if not self.isNormalised:
				sK = np.sqrt(np.pi) * sigma * sK

			#if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
			if self.isNegativeS:
				K = self.sensitivity * sK;
			else:
				K = np.sqrt(self.variance) *sK
			#if isfield(simKern2, 'isNegativeS') && (simKern2.isNegativeS == true)
			if other.isNegativeS:
				K = other.sensitivity * K
			else:
				K = np.sqrt(other.variance) * K
			return K

			
		elif isinstance(other,rbf):
			assert self.isNormalised == other.isNormalised
			assert self.inverseWidth == other.inverseWidth
			assert other.variance==1.


			dim1 = t1.shape[0]
			dim2 = t2.shape[0]
			t1 = t1 - self.delay;

			#t1Mat = t1(:, ones(1, dim2))
			t1Mat = t1*np.ones((1,dim2))
			#t2Mat = t2(:, ones(1, dim1))';
			t2Mat = np.ones((dim1,1))*t2.T

			diffT = t1Mat - t2Mat

			sigma = np.sqrt(2./self.inverseWidth)

			invSigmaDiffT = 1./sigma*diffT
			halfSigmaD_i = 0.5*sigma*self.decay

			if not self.isStationary:
				lnPart, signs = lnDiffErfs(halfSigmaD_i + t2Mat/sigma, halfSigmaD_i - invSigmaDiffT)
			else:
				lnPart, signs = lnDiffErfs(inf, halfSigmaD_i - invSigmaDiffT)

			sK = signs * np.exp(halfSigmaD_i*halfSigmaD_i - self.decay*diffT + lnPart)

			sK = 0.5 * sK

			if not self.isNormalised:
				sK = sK * np.sqrt(np.pi);
				if self.isNegativeS:
					K = sK * self.sensitivity * sigma
				else:
					K = sK * np.sqrt(self.variance) * sigma
			else:
				if self.isNegativeS:
					K = sK * self.sensitivity
				else:
					K = sK * np.sqrt(self.variance)
			return K

		else:
			raise NotImplementedError


