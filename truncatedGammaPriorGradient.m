function g = truncatedGammaPriorGradient(prior, x)

% TRUNCATEDGAMMAPRIORGRADIENT Gradient wrt x of the truncated gamma prior.

% PRIOR

if prior.lbound <= 0,
  Z = gammainc(prior.ubound / prior.b, prior.a);
else
  Z = gammainc(prior.ubound / prior.b, prior.a) - ...
      gammainc(prior.lbound / prior.b, prior.a);
end

% Compute gradient of prior
g = (x >= prior.lbound) .* (x <= prior.ubound) .* ((prior.a-1)./x-prior.b);
