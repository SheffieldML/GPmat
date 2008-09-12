function model = gpSubspaceCreate(q,d,X,y,options,dim)

model = gpCreate(q,d,X,y,options);
model.dim = dim;
model.type = 'gpSubspace';

return;