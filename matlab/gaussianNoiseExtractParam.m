function [params, names] = gaussianNoiseExtractParam(noise)

% GAUSSIANNOISEEXTRACTPARAM Extract parameters from Gaussian noise model.

% IVM

params = [noise.bias log(noise.sigma2)];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  names{i+1} = ['Log sigma^2 ' num2str(i)];
end