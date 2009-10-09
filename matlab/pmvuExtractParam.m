function [params, names] = pmvuExtractParam(model)

% PMVUEXTRACTPARAM Extract parameters from the PMVU model structure.
% FORMAT
% DESC extracts parameters from the probabilistic maximum variance unfolding
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the probabilistic maximum variance unfolding
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO pmvuCreate, pmvuExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence 2009
%
% MLTOOLS


  params = model.kappa(:)';
  fhandle = str2func([model.kappaTransform 'Transform']);
  params = fhandle(params, 'xtoa');

  counter = 0;
  for j = 1:size(model.kappa, 2)
    for i = 1:size(model.kappa, 1)
      counter = counter + 1;
      names{counter} = ['Spring ' num2str(i) ' to ' num2str(model.indices(i, j))] ;
    end
  end
end
