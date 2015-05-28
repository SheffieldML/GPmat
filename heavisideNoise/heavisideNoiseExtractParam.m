function [params, names] = heavisideNoiseExtractParam(noise)

% HEAVISIDENOISEEXTRACTPARAM Extract parameters from heaviside noise model.

% IVM

params = [invSigmoid(2*noise.eta) noise.bias];

if nargout > 1
  names{1} = 'inv logistic of eta';
  for i = 1:noise.numProcess
    names{1+i} = ['bias ' num2str(i)];
  end
end