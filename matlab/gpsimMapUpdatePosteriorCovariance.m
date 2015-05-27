function model = gpsimMapUpdatePosteriorCovariance(model)

% GPSIMMAPUPDATEPOSTERIORCOVARIANCE update the posterior covariance of f.
% FORMAT
% DESC updates the stored representation of the posterior
% covariance of f in a GPSIMMAP model. Relies on the representation
% of the kernel and the data Hessian being up to date. These can be
% forced to be up to date by using gpsimMapUpdateKernels and
% gpsimMapFunctionalUpdateW. 
% ARG model : the model for which the representation is to be
% updated.
% RETURN model : the model with the updated representation.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008
%
% SEEALSO : gpsimMapUpdateKernels, gpsimMapFunctionalUpdateW

% SHEFFIELDML


[U, Lambda] = eig(model.W);

lambda = diag(Lambda);
[lambda, order] = sort(lambda, 1, 'descend');
U = U(:, order);
% Hack to deal with non positive definite matrix.
lambda(find(lambda<0)) = 0;
ind = find(lambda)~=0;
LambdaHalf  = sparseDiag(sqrt(lambda));
LambdaHalf = LambdaHalf(ind, ind);
model.Whalf = U(:, ind)*LambdaHalf;
if size(model.Whalf, 2) == 0;
  model.Whalf = zeros(size(model.Whalf, 1), 1);
end
model.invCovf = model.invK + model.Whalf*model.Whalf';
KWhalf = model.K*model.Whalf;
inner = eye(size(model.Whalf, 2)) +model.Whalf'*KWhalf;
[innerInv, U] = pdinv(inner);
logDetInner = logdet(inner, U); 
model.covf = model.K - KWhalf*innerInv*KWhalf';
model.logDetCovf = (model.logDetK - logDetInner);                                                 
%[model.covf, U, jitter] = pdinv(model.invCovf);
%if jitter > 1e-4
%  fprintf('Warning: gpsimMapUpdatePosteriorCovariance added jitter of %2.4f\n', jitter)
%end
%model.logDetCovf = - logdet(model.invCovf, U); 
if isfield(model,'priorProtein') && ~isempty(model.priorProtein)
  nCons = length(model.priorProtein);
  consMat = zeros(size(model.covf));
  for k = 1:nCons
    ftimeIndex = find((model.priorProteinTimes(k)-model.mapt)==0);
    consMat(ftimeIndex,ftimeIndex) = model.consLambda;
  end
  hf = inv(model.invCovf+consMat);
else
  hf = model.covf;
end

model.varf = diag(hf);
