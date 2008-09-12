function model = gpSubspaceOptimise(model,varargin)

model = gpOptimise(model,varargin{:});

return;