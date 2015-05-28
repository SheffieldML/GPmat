function handle = spectrumVisualise(spectrumValues, convertFile, ...
                                    width, varargin)

% SPECTRUMVISUALISE Helper code for showing an spectrum during 2-D visualisation.
% MLTOOLS

if nargin < 3
  width = 1000;
end
cData =zeros(length(spectrumValues), width);
if nargin > 1
  spectrumValues = feval(convertFile, spectrumValues, varargin{:});
end
cData(1, 1) = -80;

handle = imagesc(cData);
