function [model, options] = modelOptimise(model, varargin)

% MODELOPTIMISE Optimise the given model.
% FORMAT
% DESC is a wrapper function that optimises a given model.
% ARG model : the model to be optimised. 
% ARG X : input value to model.
% ARG y : target value for model.
% ARG display : whether or not to display optimization values.
% ARG iters : number of iterations.
% RETURN model : the optimised model.
% RETURN options : a vector containig the number of function and gradient
% evaluations. Also the number of iterations employed.
%
% SEEALSO : modelObjective, modelGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2007
 
% MLTOOLS

if nargin < 2
  varargin = {};
end
fhandle = [model.type 'Optimise'];
if exist(fhandle)==2
  fhandle = str2func(fhandle);
  model = fhandle(model, varargin{:});
else
  if ~isfield(model, 'display')
    if length(varargin)< 3
      display = 1;
    else
      display = varargin{3};
    end    
    if length(varargin)<4
      iters = 500;
    else
      iters = varargin{4};
    end
    if length(varargin)<2 | isempty(varargin{2})
    else
      model.y = varargin{2};
    end
    if length(varargin)<1 | isempty(varargin{1})
    else
      model.X = varargin{1};
    end
  end

  options = optOptions;
  options(14) = iters;
  options(9) = 0;
  options(1) = display;
  options(2) = 1e-6;
  
  
  params = modelExtractParam(model);
  if(~isempty(params))
    if isfield(model, 'optimiser')
      optim = str2func(model.optimiser);
    else
      optim = str2func('conjgrad');
    end
    
    [params, options] = optim('modelObjective', params,  options, ...
                   'modelGradient', model);
    
    model = modelExpandParam(model, params);
  else
    warning('This Model Has No Parameters To Optimise');
  end
end


