function y = orderedNoiseOut(noise, mu, varsigma)

% ORDEREDNOISEOUT Output from ordered categorical noise model.

% NOISE

% NOISE

D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

y = zeros(size(mu));
i = noise.C-2;
index = find(mu > sum(noise.widths(1:i)));
y(index) = i+1;
for i = noise.C-3:-1:0
  index = find(mu > sum(noise.widths(1:i)) ...
               & mu <= sum(noise.widths(1:i+1)));
  y(index) = i+1;
end