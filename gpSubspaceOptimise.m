function model = gpSubspaceOptimise(model,varargin)

% GPSUBSPACEOPTIMISE
%
% COPYRIGHT : Carl Henrik Ek, 2008
  
% GPMAT

model = gpOptimise(model,varargin{:});

return;