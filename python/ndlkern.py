##SETUP
import sys
import os
import posix
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'swig', 'src'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'optimi', 'python'))
sys.path.append(os.path.join(posix.environ['HOME'], 'mlprojects', 'python', 'netlab'))
import pdb
##ENDSETUP

import math
import ndloptimi as optimi
import netlab
import numpy as np

class rbfkern():

    def __init__(self, data):
        self.type = 'rbf'
        self.inverseWidth = 1.0
        self.variance = 1.0
        self.nParams = 2
        
        # Constrains parameters positive for optimisation.
        self.transforms = optimi.transforms([0, 1], optimi.defaultConstraint('positive'))
        self.isStationary = True
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
         % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''
         
        if x2 is None:
            x = np.asmatrix(x)
            n2 = netlab.dist2(x, x)
            wi2 = (.5 * self.inverseWidth)
            k = self.variance*np.exp(-n2*wi2)
        else:
            x = np.asmatrix(x)
            x2 = np.asmatrix(x2)
            n2 = dist2(x, x2)
            wi2 = (.5 * self.inverseWidth)
            k = np.asmatrix(self.variance*np.exp(-n2*wi2))
        
        return k, n2



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
        % COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006'''

        k = np.asmatrix(np.tile(self.variance, (x.shape[0], 1)))

    def diagGradX(kern, X):
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
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009'''
        
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
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''


        g = np.zeros(1, self.nParams)
        g[0, 0] = 0.0
        g[0, 1] = covDiag.sum()
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
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009'''
        
        spacing = ''
        for i in range(numSpaces):
            spacing = spacing + ' '
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
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''

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
        % COPYRIGHT : Neil D. Lawrence, 2004--2006, 2009'''
        
        return np.array([self.inverseWidth, self.variance])

    def extractParamNames(self):
        '''% EXTRACTPARAMNAMES Extract parameter names from the RBF kernel structure.
        % FORMAT
        % DESC Extract parameter names from the radial basis
        % function kernel structure.
        % ARG kern : the kernel structure containing the parameters to be
        % extracted.
        % RETURN names : cell array of strings giving names to the parameters.'''
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
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''
        
        gX = np.empty((X2.shape[0], X2.shape[1], X.shape[0]))
        for i in range(X.shape[0]):
            gX[:, :, i] = self.gradXpoint(X[i, :], X2)
        return gX


    def gradXpoint(self, x, X2):
        '''% GRADXPOINT Gradient with respect to one point of x.'''

        gX = np.zeros(X2.shape)
        n2 = netlab.dist2(X2, x)
        wi2 = (.5 * self.inverseWidth);
        rbfPart = self.variance*exp(-n2*wi2);
        for i in range(x.shape[1]):
            gX[:, i] = self.inverseWidth*(X2[:, i] - x[i])*rbfPart
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
        % COPYRIGHT : Neil D. Lawrence, 2004-2006, 2009'''

        if X2 is None:
            k, dist2xx = self.compute(X)
        else:
            k, dist2xx = compute(X, X2)
        g = np.zeros((1, self.nParams))
        g[0, 0] = - .5*(np.asarray(covGrad)*np.asarray(k)*np.asarray(dist2xx)).sum()
        g[0, 1] =  (np.asarray(covGrad)*np.asarray(k)).sum()/self.variance;
        return g
