#!/usr/bin/env python

# For creating files for new model. Run as 

# modelGenerator.py {model short name} {model long name} {toolbox name} [CopyrightName YearOne YearTwo etc]

# A detailed model description can be given in {model short name}.txt
import os
import sys
import re
import shutil
import time
import string


if len(sys.argv) < 4:
    raise "There should be three input arguments"
prefix = sys.argv[1]
longName = sys.argv[2]
toolboxName = sys.argv[3]
toolboxCapName = toolboxName.upper()
ucPrefix = prefix.upper()
year = time.strftime('%Y')

if len(sys.argv)>4:
    # copyright info provided
    copyRightText = sys.argv[4]
    if len(sys.argv)>5:
        for argNo in range(len(sys.argv)-5):
            copyRightText += ', ' + sys.argv[argNo+5]
    else:
        copyRightText + ', ' + year
    copyRightText = '\n%\n% COPYRIGHT : ' + copyRightText
else:
    copyRightText = ''
    
file = prefix + '.txt'
if os.path.exists(file):
    fileHandle = open(file, 'r')
    modelDescription = fileHandle.read()
    fileHandle.close()
    modelDescription = '\n' + modelDescription.strip() + '\n%'
else:
    modelDescription = ''

# Dictionary of the files and their contents.
files = {prefix + 'Out.m' :
 '''function y = ''' + prefix + '''Out(model, x)

% ''' + ucPrefix + '''OUT Compute the output of a ''' + ucPrefix + ''' model given the structure and input X.
% FORMAT
% DESC computes the model parameters for the ''' + longName + '''
% model given inputs associated with rows and columns.
% ARG model : the model structure for which the output is computed.
% ARG x : the input data.
% RETURN y : the output results.
%
% SEEALSO : ''' + prefix + '''Create, modelCompute, modelCreate, ''' + prefix + '''ExpandParam, ''' + prefix + '''ExtractParam''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'Deconstruct.m' :
'''function ''' + prefix + '''Info = ''' + prefix + '''Deconstruct(model)

% ''' + ucPrefix + '''DECONSTRUCT break ''' + ucPrefix + ''' in pieces for saving.
% FORMAT
% DESC takes an ''' + longName + ''' model structure and breaks it into component
% parts for saving. 
% ARG model : the model that needs to be saved.
% RETURN ''' + prefix + '''Info : a structure containing the other information
% from the ''' + longName + ''': what the sparse approximation is, what the inducing
% variables are.
%
% SEEALSO : ''' + prefix + '''Create, ''' + prefix + '''Reconstruct''' + copyRightText + '''
 
% ''' + toolboxCapName + '''

''',
prefix + 'Reconstruct.m' :
'''function model = ''' + prefix + '''Reconstruct(''' + prefix + '''Info, y)

% ''' + ucPrefix + '''RECONSTRUCT Reconstruct an ''' + longName + ''' from component parts.
% FORMAT
% DESC takes component parts of an ''' + longName + ''' model and reconstructs the
% ''' + longName + ''' model. The component parts are normally retrieved from a
% saved file.
% ARG ''' + prefix + '''Info : the active set and other information stored in a structure.
% ARG y : the output target training data for the ''' + longName + '''.
% RETURN model : an ''' + longName + ''' model structure that combines the component
% parts.
%
% SEEALSO : ''' + prefix + '''Create, ''' + prefix + '''Reconstruct''' + copyRightText + '''
 
% ''' + toolboxCapName + '''

''',
prefix + 'LogLikelihood.m' :
'''function ll = ''' + prefix + '''LogLikelihood(model)

% ''' + ucPrefix + '''LOGLIKELIHOOD Log likelihood of ''' + ucPrefix + ''' model.
% FORMAT
% DESC computes the log likelihood of  the ''' + longName + ''' model.
% ARG model : the model structure for which log likelihood is being computed.
% RETURN ll : the computed log likelihood.
%
% SEEALSO : ''' + prefix + '''Create, ''' + prefix + '''LogLikeGradients, modelLogLikelihood''' + copyRightText + '''
 
% ''' + toolboxCapName + '''

''',
prefix + 'Display.m' : 
'''function ''' + prefix + '''Display(model, spacing)

% ''' + ucPrefix + '''DISPLAY Display parameters of the ''' + ucPrefix + ''' model.
% FORMAT
% DESC displays the parameters of the ''' + longName + '''
% model and the model type to the console.
% ARG model : the model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG model : the model to display.
% ARG spacing : how many spaces to indent the display of the model by.
%
% SEEALSO : ''' + prefix + '''Create, modelDisplay''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'ExpandParam.m' :
'''function model = ''' + prefix + '''ExpandParam(model, params)

% ''' + ucPrefix + '''EXPANDPARAM Create model structure from ''' + ucPrefix + ''' model\'s parameters.
% FORMAT
% DESC returns a ''' + longName + ''' model structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG model : the model structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% model structure.
% RETURN model : model structure with the given parameters in the
% relevant locations.
%
% SEEALSO : ''' + prefix + '''Create, ''' + prefix + '''ExtractParam, modelExpandParam''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'ExtractParam.m' :
'''function [params, names] = ''' + prefix + '''ExtractParam(model)

% ''' + ucPrefix + '''EXTRACTPARAM Extract parameters from the ''' + ucPrefix + ''' model structure.
% FORMAT
% DESC extracts parameters from the ''' + longName + '''
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the ''' + longName + '''
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO ''' + prefix + '''Create, ''' + prefix + '''ExpandParam, modelExtractParam, scg, conjgrad''' + copyRightText + '''
%
% ''' + toolboxCapName + '''

''',
prefix + 'LogLikeGradients.m' : 
'''function g = ''' + prefix + '''LogLikeGradients(model)

% ''' + ucPrefix + '''LOGLIKEGRADIENTS Gradient of ''' + ucPrefix + ''' model log likelihood with respect to parameters.
% FORMAT
% DESC computes the gradient of the ''' + longName + '''
% model\'s log likelihood with respect to the parameters.
% ARG model : model structure for which gradients are being
% computed.
% RETURN g : the returned gradients. 
%
% SEEALSO ''' + prefix + '''Create, ''' + prefix + '''LogLikelihood, modelLogLikeGradients ''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'ParamInit.m' : 
'''function model = ''' + prefix + '''ParamInit(model)

% ''' + ucPrefix + '''PARAMINIT ''' + ucPrefix + ''' model parameter initialisation.
% FORMAT
% DESC initialises the ''' + longName + '''
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : ''' + prefix + '''Create, modelCreate, modelParamInit''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'Create.m' : 
'''function model = ''' + prefix + '''Create(inputDim, outputDim, options)

% ''' + ucPrefix + '''CREATE Create a ''' + ucPrefix + ''' model.''' + modelDescription + '''
% FORMAT
% DESC creates a ''' + longName + '''
% model structure given an options structure. 
% ARG inputDim : the input dimension of the model.
% ARG outputDim : the output dimension of the model.
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : ''' + prefix + '''Options, ''' + prefix + '''ParamInit, modelCreate''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'Options.m' :
'''function options = ''' + prefix + '''Options

% ''' + ucPrefix + '''OPTIONS Create a default options structure for the ''' + ucPrefix + ''' model.
% FORMAT
% DESC creates a default options structure for the ''' + longName + ''' model
% structure.
% RETURN options : the default options structure.
%
% SEEALSO : ''' + prefix + '''Create, modelOptions''' + copyRightText + '''

% ''' + toolboxCapName + '''

''',
prefix + 'OutputGradX.m' :
''' function g = ''' + prefix + '''OutputGradX(model, X)

% ''' + ucPrefix + '''OUTPUTGRADX Evaluate derivatives of a ''' + ucPrefix + ''' model\'s output with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an ''' + longName + ''' model with
% respect to the inputs to the model. 
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs.
%
% SEEALSO : ''' + prefix + '''OutputGrad, modelOutputGradX''' + copyRightText + '''

% ''' + toolboxCapName + '''
''',
prefix + 'OutputGrad.m' :
''' function g = ''' + prefix + '''OutputGrad(model, X)

% ''' + ucPrefix + '''OUTPUTGRAD Evaluate derivatives of ''' + ucPrefix + ''' model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of a ''' + longName + ''' model
% outputs with respect to the parameters of the ''' + longName + '''
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the ''' + longName + '''
% with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : ''' + prefix + '''Create, ''' + prefix + '''LogLikeGradients''' + copyRightText + '''

% ''' + toolboxCapName + '''

'''}

for file in files.keys():
    if os.path.exists(file):
        foundToolboxName = 0
        justFoundToolboxName = 0
        newlines = files[file].split('\n')
        funcLine = newlines[0] + '\n'
        writeText = ''
        fileHandle = open(file, 'r')
        lines = fileHandle.readlines()
        fileHandle.close()
        for line in lines:
            if re.findall(re.compile(r'^%\s*'+toolboxCapName+'\s*'), line):
                foundToolboxName = 1
                justFoundToolboxName = 1
            elif justFoundToolboxName:
                if not re.findall(re.compile(r'^\s*$'), line):
                    writeText += line
                    justFoundToolboxName = 0
            elif foundToolboxName:
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
