function model = gpSubspaceCreate(q,d,X,y,options,dim)

% GPSUBSPACECREATE 
%
% COPYRIGHT : Carl Henrik Ek, 2008

% GPMAT
  
model = gpCreate(q,d,X,y,options);
model.dim = dim;
model.type = 'gpSubspace';

return;