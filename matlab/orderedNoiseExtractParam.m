function [params, names] = orderedNoiseExtractParam(noise)

% ORDEREDNOISEEXTRACTPARAM Extract parameters from ordered categorical noise model.

% IVM
params = [noise.bias noise.widths(:)'];
if nargout > 1
  for i = 1:noise.numProcess
    names{i} = ['bias ' num2str(i)];
  end
  for i = 1:noise.C-2
    names{noise.numProcess+i} = ['width ' num2str(i)];
  end
end