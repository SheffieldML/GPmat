function model = gpSubspaceOptimise(model,varargin)

% GPSUBSPACEOPTIMISE
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek


model = gpOptimise(model,varargin{:});

return;