function model = gpSubspaceCreate(q,d,X,y,options,dim)

% GPSUBSPACECREATE
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek

  
model = gpCreate(q,d,X,y,options);
model.dim = dim;
model.type = 'gpSubspace';

return;