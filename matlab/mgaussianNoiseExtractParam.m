function [params, names] = mgaussianNoiseExtractParam(noise)

% MGAUSSIANNOISEEXTRACTPARAM Extract parameters from Variable variance Gaussian noise model.

% IVM


params = [noise.bias noise.sigma2];


if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  for i = noise.numProcess+1:2*noise.numProcess
    names{i} = ['sigma^2 ' num2str(i)];
  end
end