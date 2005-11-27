function handle = imageVisualise(imageVals, imageSize, transpose, negative, ...
				 scale)

% IMAGEVISUALISE Helper code for showing an image during 2-D visualisation.

% MLTOOLS

if nargin < 3
  transpose = 1;
end
if nargin< 4
  negative = 0;
end
if nargin < 5
  scale = 1;
end
if negative
  imageVals = -imageVals;
end
imageData = reshape(imageVals, imageSize(1), imageSize(2));
if transpose
  imageData = imageData';
end
if scale
  handle = imagesc(imageData);
else
  handle = image(imageData);
end