function model = gpSubspaceOptimise(model,varargin)

% GPSUBSPACEOPTIMISE
  
% GP

model = gpOptimise(model,varargin{:});

return;