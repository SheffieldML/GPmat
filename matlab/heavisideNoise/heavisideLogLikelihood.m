function L = heavisideLogLikelihood(noise, mu, varsigma, y)

% HEAVISIDELOGLIKELIHOOD Log-likelihood of data under heaviside noise model.

% IVM

% limVal = 36;
% fact = sqrt(2)/2;
% D = size(y, 2);
% L = 0;
% for i = 1:D
%   u = mu(:, i) + noise.bias(i);
%   u = (y(:, i).*u)./(sqrt(varsigma(:, i)));
%   index = find(u > limVal);
%   L = L + length(index)*log(noise.eta);
%   u(index) = [];
%   u2 = u.*u;   
%   L = L - sum(.5*u2) + sum(log((1-2*noise.eta)...
%                                *.5*erfcx(-fact*u)...
%                                +exp(.5*u2)*noise.eta));
% end
L = sum(sum(log(heavisideLikelihood(noise, mu, varsigma, y))));  