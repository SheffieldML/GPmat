function model = modelOptimise(model, varargin)

% MODELOPTIMISE Optimise the given model.

% MLTOOLS

if nargin < 2
  varargin = {}
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
end

options = optOptions;
options(14) = iters;
options(9) = 1;
options(1) = display;



params = modelExtractParam(model);
if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

params = optim('modelObjective', params,  options, ...
               'modelGradient', model);

model = modelExpandParam(model, params);


