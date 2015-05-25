function model = mogProject(model, dimension)
 
% MOGPROJECT Project a mixture of Gaussians to a low dimensional space.
% FORMAT
% DESC projects a mixture of Gaussians down to a lower dimensional space
% (typically for visualisation).
% ARG model : the mixture of Gaussians to project down.
% ARG dimension : the dimension to project to.
% RETURN model : the reduced dimensional mixture of Gaussians.
% 
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% SEEALSO : mogCreate, mogTwoDPlot

% MLTOOLS
  
  
[m, C] = mogMeanCov(model);
[U, V] = eig(C);

v = diag(V);
[v, ind] = sort(v);
ind = ind(end:-1:1);
U = U(:, ind);

Uq = U(:, 1:dimension);

model.Y = model.Y*Uq;
model.mean = model.mean*Uq;
model.d = dimension;

switch model.covtype
 case 'ppca'
  for i = 1:model.m
    model.W{i} = (model.W{i}'*Uq)';
    model.U{i} = sqrt(model.sigma2(i))*eye(model.d);
    for j = 1:model.q
      model.U{i} = cholupdate(model.U{i}, model.W{i}(:, j));
    end
  end

end
