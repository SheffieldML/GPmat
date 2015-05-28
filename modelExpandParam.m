function model = modelExpandParam(model, params, dim)

% MODELEXPANDPARAM Update a model structure with parameters.
% FORMAT
% DESC returns a model structure filled with the parameters in the
% given vector. This is used as a helper function to enable
% parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG model : the model structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% model structure.
% RETURN model : model structure with the given parameters in the
% relevant locations.
%
% SEEALSO : modelExtractParam, scg, conjgrad
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Cark Henrik Ek, 2007

% MLTOOLS

if isfield(model, 'paramGroups')
  params = params*model.paramGroups';
end

fhandle = str2func([model.type 'ExpandParam']);
if(nargin<3)
  model = fhandle(model, params);
else
  model = fhandle(model,params,dim);
end
