function [params, names] = noiseExtractParam(noise)

% NOISEEXTRACTPARAM Extract the noise model's parameters.

% NOISE

fhandle = str2func([noise.type 'NoiseExtractParam']);
if nargout < 2
  params = fhandle(noise);
else
  [params, names] = fhandle(noise);
end


% Check if parameters are being optimised in a transformed space.
if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    fhandle = str2func([noise.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end