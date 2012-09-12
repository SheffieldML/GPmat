function [mean, vars] = gpnddisimSamplePredictions(model, paramsamples, t_pred, N_samples),

% t_pred = (((0:100)/100*sqrt(1280)).^2 + 300)';

if nargin < 4,
  N_samples = 10;
end

r = zeros(size(paramsamples, 1)*N_samples, 2*length(t_pred));
baseind = 0;
for k=1:size(paramsamples, 1),
  m = modelExpandParam(model, paramsamples(k, :));
  [~, mu, C] = gpnddisimPredict(m, t_pred, 1, 0);
  C = 0.5*(C+C');
  [V,D]=eig(C);
  if min(diag(D)) < -1e-8,
    fprintf('min(D) = %g\n', min(diag(D)));
  end
  r(baseind+1:baseind+N_samples, :) = mvnrnd(mu', C-1.1*min(0, min(diag(D))) * eye(size(C)), N_samples);
  baseind = baseind + N_samples;
end
mean = r;
vars = [];
