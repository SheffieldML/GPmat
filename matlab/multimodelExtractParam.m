function [params, names] = multimodelExtractParam(model)

% MULTIMODELEXTRACTPARAM Extract parameters from the MULTIMODEL model structure.
% FORMAT
% DESC extracts parameters from the multi-task learning wrapper
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the multi-task learning wrapper
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO multimodelCreate, multimodelExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% MLTOOLS

[params, names] = modelExtractParam(model.comp{1});
