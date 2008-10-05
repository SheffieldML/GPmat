function model = gpSubspaceCreate(q,d,X,y,options,dim)

% GPSUBSPACECREATE 
  
% GP
  
model = gpCreate(q,d,X,y,options);
model.dim = dim;
model.type = 'gpSubspace';

return;