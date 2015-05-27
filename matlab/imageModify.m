function handle = imageModify(handle, imageValues, imageSize, transpose, negative, ...
		     scale)

% IMAGEMODIFY Helper code for visualisation of image data.
% FORMAT
% DESC is a helper function for visualising image data using latent
% variable models.
% ARG handle : the handle of the image data.
% ARG imageValues : the values to set the image data to.
% ARG imageSize : the size of the image.
% ARG transpose : whether the resized image needs to be transposed
% (default 1, which is yes).
% ARG negative : whether to display the negative of the image
% (default 0, which is no).
% ARG scale : dummy input, to maintain compatability with
% IMAGEVISUALISE.
% RETURN handle : a the handle to the image data.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
%
% SEEALSO : imageVisualise, fgplvmResultsDynamic

% MLTOOLS


if nargin < 4
  transpose = 1;
end
if nargin< 5
  negative = 0;
end
if negative
  imageValues = -imageValues;
end
if transpose
  set(handle, 'CData', reshape(imageValues(1:imageSize(1)*imageSize(2)), imageSize(1), imageSize(2))');
else
  set(handle, 'CData', reshape(imageValues(1:imageSize(1)*imageSize(2)), imageSize(1), imageSize(2)));
end
