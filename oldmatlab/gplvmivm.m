function [activeSet, Kstore] = gplvmivm(X, theta, d)

% GPLVMIVM THis code implements active set selection (via the IVM) for the GPLVM

numData = size(X, 1);
% Create storage for the kernel
Kstore = zeros(numData, d);
% Initialise Indices
J = 1:numData;
activeSet = [];

% Note the last parameter theta(end) is the precision.

% Initialise parameters
diagK = kerneldiag(X, theta);
diagA = diagK;
%h = zeros(numData, 1);

[infoChange, indexSelect] = max(diagK);
i = J(indexSelect);

k = 1;

% Update L, M, h and A with selected point
L = sqrt(1 + theta(end)*diagK(i));
Kstore(:, 1) = kernel(X, X(i, :), theta);
Kstore(i, 1) = Kstore(i, 1) + 1/theta(end);

M = (1/L*(sqrt(theta(end))*Kstore(:, 1)))';
diagA = diagA - M'.*M';
%h = h + alpha(i)*L*1/sqrt(theta(end))*M';

% Remove point from the non-active set and place in the active.
J(indexSelect) = [];
activeSet = [activeSet; i];

% Update nu
infoChange(1) = NaN;

% Record Likelihood
%logLikelihoods = log(cummGaussian(z));
%logLikelihood = sum(logLikelihoods);

for k = 2:d
  delta = -.5*log2(1-diagA(J));
  [infoChange(k), indexSelect] = max(delta);
  i = J(indexSelect);
  
  % Update L, M, h and A with selected point
  lvec = sqrt(theta(end))*M(:, i);
  l = sqrt(1 + theta(end)*diagK(i) - lvec'*lvec);
  Kstore(:, k) = kernel(X, X(i, :), theta);
  Kstore(i, k) = Kstore(i, k) + 1/theta(end); % introduce diagonal term
  muvec = (1/l)*sqrt(theta(end))*(Kstore(:, k) - M'*M(:, i)); % Neil's
  
  L = [L zeros(k-1, 1); lvec' l];
  M = [M; muvec'];
  diagA = diagA - muvec.*muvec;
  %h = h + alpha(i)*l*1/sqrt(theta(end))*muvec;
  
  
  % Remove point from the non-active set and place in the active.
  J(indexSelect) = [];
  activeSet = [activeSet; i];
  
  
%  logLikelihoods = log(cummGaussian(z));
%  dLogLikelihood(k) = sum(logLikelihoods);
%  falsePositives(k) = sum(sign(h)~=y & y==1);
%  trueNegatives(k) = sum(sign(h)~=y & y==-1);
end






