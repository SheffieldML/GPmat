function g = mixturePriorGradient(prior, x)

% MIXTUREPRIORGRADIENT Gradient wrt x of the mixture prior.

% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

grads = cellfun(@(p) priorGradient(p, x), prior.comp, 'UniformOutput', false);
logprobs = cellfun(@(p) priorLogProb(p, x), prior.comp) + log(prior.weights);

normps = exp(logprobs - max(logprobs));
g = zeros(size(x));
for k=1:length(normps),
  g = g + normps(k) * grads{k};
end
g = g / sum(normps);
