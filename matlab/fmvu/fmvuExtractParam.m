function [params, names] = fmvuExtractParam(model)

% FMVUEXTRACTPARAM Extract parameters from the FMVU model structure.
% FORMAT
% DESC extracts parameters from the fast maximum variance unfolding
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the fast maximum variance unfolding
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO fmvuCreate, fmvuExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS


if nargout > 1
  returnNames = true;
else
  returnNames = false;
end  

params = [model.X(:)' log(model.kappa(:))'];
if returnNames
  for i = 1:model.N
    for j = 1:model.q
      Xnames{i, j} = ['X(' num2str(i) ', ' num2str(j) ')'];
    end
    for j = 1:model.k
      Kapnames{i, j} = ['kappa(' num2str(i) ', ' num2str(j) ')'];
    end
  end
  names = {Xnames{:}, Kapnames{:}};
end
