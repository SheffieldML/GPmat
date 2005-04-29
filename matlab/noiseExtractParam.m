function [params, names] = noiseExtractParam(noise)

% NOISEEXTRACTPARAM Extract the noise model's parameters.

% NOISE

if nargout < 2
  params = feval([noise.type 'NoiseExtractParam'], noise);
else
  [params, names] = feval([noise.type 'NoiseExtractParam'], noise);
end


% Check if parameters are being optimised in a transformed space.
if isfield(noise, 'transforms')
  for i = 1:length(noise.transforms)
    index = noise.transforms(i).index;
    params(index) = feval([noise.transforms(i).type 'Transform'], ...
              params(index), 'xtoa');
  end
end