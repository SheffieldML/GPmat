function Y = expvarMeanOut(model, X);

% EXPVARMEANOUT Output of an EXPVARMEAN model.
% FORMAT 
% DESC gives the output of the mean of a variational approximation
% to a given kernel. This is for use with the exponentiated kernel
% when making a variational approximation to an exponentiated
% Gaussian process.
% ARG model : the model for which the output is required.
% ARG X : the input data for which the output is required.
% RETURN Y : the output.
%
%
% SEEALSO :  modelOut, expvarMeanCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

Y = exp(0.5*kernDiagCompute(model.kern.argument, X));
