function l = truncatedGammaPriorLogProb(prior, x)

% TRUNCATEDGAMMAPRIORLOGPROB Log probability of truncated gamma prior.

% PRIOR

% Compute log prior

% First compute the normalising factor due to truncation
if prior.lbound <= 0,
  Z = gammainc(prior.ubound * prior.b, prior.a);
else
  Z = gammainc(prior.ubound * prior.b, prior.a) - ...
      gammainc(prior.lbound * prior.b, prior.a);
end
  
D = length(x);
l = sum(log((x >= prior.lbound) .* (x <= prior.ubound) * ...
            1 / Z)) + ...
    D*prior.a*log(prior.b)-D*gammaln(prior.a)+(prior.a-1)*sum(log(x))-prior.b*sum(x);
