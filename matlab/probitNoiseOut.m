function y = probitNoiseOut(noise, mu, varsigma)

% PROBITNOISEOUT Output from probit noise model.

% IVM
D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = mu > 0;
