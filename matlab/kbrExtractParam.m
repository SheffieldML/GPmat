function params = kbrExtractParam(model);

% KBREXTRACTPARAM Extract parameters from the KBR model structure.
% FORMAT
% DESC extracts parameters from the kernel based regression
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the kernel based regression
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO kbrCreate, kbrExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MLTOOLS


params = [model.A(:)' model.bias];