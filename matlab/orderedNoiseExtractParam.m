function [params, names] = orderedNoiseExtractParam(noise)

% ORDEREDNOISEEXTRACTPARAM Extract parameters from ordered categorical noise model.

% IVM

params = [invSigmoid(noise.C*noise.eta) noise.bias log(noise.widths(:))'];

if nargout > 1
  names{1} = 'inv logistic of eta';
  for i = 1:noise.numProcess
    names{1+i} = ['bias ' num2str(i)];
  end
  for i = 1:noise.C-2
    names{noise.numProcess+i+1} = ['log width ' num2str(i)];
  end
end