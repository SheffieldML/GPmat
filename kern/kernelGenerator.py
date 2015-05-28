#!/usr/bin/env python

# For creating files for new kernels. Run as 

# kernelGenerator.py {kernel short name} {kernel long name} [CopyrightName YearOne YearTwo etc]

# A detailed kernel description can be given in {kernel short name}.txt
import os
import sys
import re
import shutil
import time
import string


if len(sys.argv) < 3:
    raise "There should be two input arguments"
prefix = sys.argv[1]
longName = sys.argv[2]
ucPrefix = prefix.upper()
year = time.strftime('%Y')

if len(sys.argv)>3:
    # copyright info provided
    copyRightText = sys.argv[3]
    if len(sys.argv)>4:
        for argNo in range(len(sys.argv)-4):
            copyRightText += ', ' + sys.argv[argNo+4]
    else:
        copyRightText + ', ' + year
    copyRightText = '\n%\n% COPYRIGHT : ' + copyRightText
else:
    copyRightText = ''
    
file = prefix + '.txt'
if os.path.exists(file):
    fileHandle = open(file, 'r')
    kernDescription = fileHandle.read()
    fileHandle.close()
    kernDescription = '\n' + kernDescription.strip() + '\n%'
else:
    kernDescription = ''

# Dictionary of the files and their contents.
files = {prefix + 'KernCompute.m' :
 '''function K = ''' + prefix + '''KernCompute(kern, x, x2)

% ''' + ucPrefix + '''KERNCOMPUTE Compute the ''' + ucPrefix + ''' kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the ''' + longName + '''
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the ''' + longName + '''
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : ''' + prefix + '''KernParamInit, kernCompute, kernCreate, ''' + prefix + '''KernDiagCompute''' + copyRightText + '''

% KERN

''',
prefix + 'KernDiagCompute.m' : 
'''function k = ''' + prefix + '''KernDiagCompute(kern, x)

% ''' + ucPrefix + '''KERNDIAGCOMPUTE Compute diagonal of ''' + ucPrefix + ''' kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the ''' + longName + ''' kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : ''' + prefix + '''KernParamInit, kernDiagCompute, kernCreate, ''' + prefix + '''KernCompute''' + copyRightText + '''

% KERN

''',
prefix + 'KernDiagGradX.m' :
'''function gX = ''' + prefix + '''KernDiagGradX(kern, X)

% ''' + ucPrefix + '''KERNDIAGGRADX Gradient of ''' + ucPrefix + ''' kernel\'s diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the ''' + longName + ''' kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : ''' + prefix + '''KernParamInit, kernDiagGradX, ''' + prefix + '''kernGradX''' + copyRightText + '''

% KERN

''',
prefix + 'KernDiagGradient.template' :
''' function g = ''' + prefix + '''KernDiagGradient(kern, x, covDiag)

% ''' + ucPrefix + '''KERNDIAGGRADIENT Compute the gradient of the ''' + ucPrefix + ''' kernel\'s diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% ''' + longName + ''' kernel matrix with respect to the parameters of the kernel. The
% parameters\' gradients are returned in the order given by the
% ''' + prefix + '''KernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% ''' + prefix + '''KernExtractParam.
%
% SEEALSO : ''' + prefix + '''KernParamInit, kernDiagGradient, ''' + prefix + '''KernExtractParam, ''' + prefix + '''KernGradient''' + copyRightText + '''

% KERN

''',
prefix + 'KernDisplay.m' : 
'''function ''' + prefix + '''KernDisplay(kern, spacing)

% ''' + ucPrefix + '''KERNDISPLAY Display parameters of the ''' + ucPrefix + ''' kernel.
% FORMAT
% DESC displays the parameters of the ''' + longName + '''
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : ''' + prefix + '''KernParamInit, modelDisplay, kernDisplay''' + copyRightText + '''

% KERN

''',
prefix + 'KernExpandParam.m' :
'''function kern = ''' + prefix + '''KernExpandParam(kern, params)

% ''' + ucPrefix + '''KERNEXPANDPARAM Create kernel structure from ''' + ucPrefix + ''' kernel\'s parameters.
% FORMAT
% DESC returns a ''' + longName + ''' kernel structure filled with the
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
% SEEALSO : ''' + prefix + '''KernParamInit, ''' + prefix + '''KernExtractParam, kernExpandParam''' + copyRightText + '''

% KERN

''',
prefix + 'KernExtractParam.m' :
'''function [params, names] = ''' + prefix + '''KernExtractParam(kern)

% ''' + ucPrefix + '''KERNEXTRACTPARAM Extract parameters from the ''' + ucPrefix + ''' kernel structure.
% FORMAT
% DESC extracts parameters from the ''' + longName + '''
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field \'transforms\' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the ''' + longName + '''
% kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field \'transforms\' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO ''' + prefix + '''KernParamInit, ''' + prefix + '''KernExpandParam, kernExtractParam, scg, conjgrad''' + copyRightText + '''
%
% KERN

''',
prefix + 'KernGradX.m' : 
'''function gX = ''' + prefix + '''KernGradX(kern, x, x2)

% ''' + ucPrefix + '''KERNGRADX Gradient of ''' + ucPrefix + ''' kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the ''' + longName + '''
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
% DESC computes the gradident of the ''' + longName + '''
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
% SEEALSO ''' + prefix + '''KernParamInit, kernGradX, ''' + prefix + '''KernDiagGradX''' + copyRightText + '''

% KERN

''',
prefix + 'KernGradient.m' :
'''function g = ''' + prefix + '''KernGradient(kern, x, x2, covGrad)

% ''' + ucPrefix + '''KERNGRADIENT Gradient of ''' + ucPrefix + ''' kernel\'s parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% '''+ longName + '''
% kernel\'s parameters. As well as the kernel structure and the
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
% SEEALSO ''' + prefix + '''KernParamInit, kernGradient, ''' + prefix + '''KernDiagGradient, kernGradX''' + copyRightText + '''

% KERN

''',
prefix + 'KernParamInit.m' : 
'''function kern = ''' + prefix + '''KernParamInit(kern)

% ''' + ucPrefix + '''KERNPARAMINIT ''' + ucPrefix + ''' kernel parameter initialisation.''' + kernDescription + '''
% FORMAT
% DESC initialises the ''' + longName + '''
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit''' + copyRightText + '''

% KERN

'''};


for file in files.keys():
    if os.path.exists(file):
        foundKern = 0
        justFoundKern = 0
        newlines = files[file].split('\n')
        funcLine = newlines[0] + '\n'
        writeText = ''
        fileHandle = open(file, 'r')
        lines = fileHandle.readlines()
        fileHandle.close()
        for line in lines:
            if re.findall(re.compile(r'^%\s*KERN\s*'), line):
                foundKern = 1
                justFoundKern = 1
            elif justFoundKern:
                if not re.findall(re.compile(r'^\s*$'), line):
                    writeText += line
                    justFoundKern = 0
            elif foundKern:
                writeText += line
            elif re.findall(re.compile(r'^\s*function.*'), line):
                funcLine = line
        for line in newlines:
            if re.findall(re.compile(r'^\s*function.*'), line):
                commentText = funcLine + '\n'
            else:
                commentText += line + '\n'
        writeText = commentText + writeText
        fileHandle = open(file, 'w')
        fileHandle.write(writeText)
        fileHandle.close()
        
    else:
        fileHandle = open(file, 'w')
        fileHandle.write(files[file])
        fileHandle.close()
