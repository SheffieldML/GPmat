function model = expvarMeanCreate(inputDim, outputDim, options)

% EXPVARMEANCREATE Creates an the mean function for the EXP kernel.
% FORMAT
% DESC creates a model for returning the first moment of an
% 'exponentiated Gaussian process'. If the output of a Gaussian
% process is exponentiated the resulting process is no longer
% Gaussian. However, it can be approximated by a Gaussian process
% by matching its first and second moments. This function returns
% the mean of that approximating process for a given kernel
% (specified in the options vector). It should be used in tandem
% with the EXP kernel for approximating these Gaussian processes.
% ARG inputDimension : dimension of input to function.
% ARG outputDim : dimension of output from mean function data.
% ARG options : options structure. The structure contains the type
% of kernel that the function is based on. A set of default options 
% are given by the file expvarMeanOptions.
% RETURN model : model structure containing the mapping.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : expvarMeanOptions

% SHEFFIELDML

model.type = 'expvarMean';
if isstruct(options.kern) 
  model.kern = options.kern;
else
  model.kern = kernCreate(inputDim, options.kern);
end
model.q = inputDim;
model.d = outputDim;
model.numParams = model.kern.nParams;
