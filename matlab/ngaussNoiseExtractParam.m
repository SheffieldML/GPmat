function [params, names] = ngaussNoiseExtractParam(noise)

% NGAUSSNOISEEXTRACTPARAM Extract parameters from noiseless Gaussian noise model.

% NOISE

% NOISE


params = [noise.bias];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
end