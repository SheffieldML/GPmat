function [params, names] = probitNoiseExtractParam(noise)

% PROBITNOISEEXTRACTPARAM Extract parameters from probit noise model.

% NOISE

params = [noise.bias];

if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
end