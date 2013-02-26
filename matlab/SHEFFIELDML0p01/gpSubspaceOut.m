function y = gpSubspaceOut(model,x)

% GPSUBSPACEOUT
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek

  
y = NaN.*ones(size(x,1),length(model.dim));
y(:,find(model.dim)) = gpOut(model,x);

return;