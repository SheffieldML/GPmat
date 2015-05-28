function handle = spectrumModify(handle, spectrumValues, convertFile, width,varargin);

% SPECTRUMMODIFY Helper code for visualisation of spectrum data.

% MLTOOLS

if nargin < 4
  width = 500;
end
persistent counter
persistent spectrum
if isempty(counter)
  counter = 1;
  spectrum = zeros(size(spectrumValues, 2), width);
end
if counter > width-1
  counter = 0;
  for i = 1:size(spectrum, 2)
    a = ep2wave(spectrum(:, i), 0.02, 'randn');
    if i == 1
      signa = a;
    else
      signa = [signa; zeros(floor(size(a, 1)/2), 1)];
      signa(end-size(a, 1)+1:end)=signa(end-size(a, 1)+1:end)+a;
    end
  end
  disp('I''m free')
  signa = signa/max(abs(signa));
  wavplay(signa, 22050);
end
counter = counter + 1;
x = repmat(counter,size(spectrumValues));
y = 1:length(spectrumValues);
%spectrumValues = diffrep2Ratemap(spectrumValues)';
spectrumValues = spectrumValues';
spectrum(:, counter) = spectrumValues;
%a = ep2wave(spectrumValues,0.02, 'randn');
%a = a/max(abs(a));
%wavplay(a, 22050);
if nargin > 2
  spectrumValues = feval(convertFile, spectrumValues, varargin{:});
end
cData = get(handle, 'CData');
cData(:, counter) = spectrumValues;
set(handle, 'CData', cData);
%disp(spectrumValues)
