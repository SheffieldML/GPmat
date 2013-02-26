function handle = imageModify(handle, imageValues, imageSize, transpose, negative, ...
		     scale)

% IMAGEMODIFY Helper code for visualisation of image data.
%
%	Description:
%
%	HANDLE = IMAGEMODIFY(HANDLE, IMAGEVALUES, IMAGESIZE, TRANSPOSE,
%	NEGATIVE, SCALE) is a helper function for visualising image data
%	using latent variable models.
%	 Returns:
%	  HANDLE - a the handle to the image data.
%	 Arguments:
%	  HANDLE - the handle of the image data.
%	  IMAGEVALUES - the values to set the image data to.
%	  IMAGESIZE - the size of the image.
%	  TRANSPOSE - whether the resized image needs to be transposed
%	   (default 1, which is yes).
%	  NEGATIVE - whether to display the negative of the image (default
%	   0, which is no).
%	  SCALE - dummy input, to maintain compatability with
%	   IMAGEVISUALISE.
%	
%
%	See also
%	IMAGEVISUALISE, FGPLVMRESULTSDYNAMIC


%	Copyright (c) 2003, 2004, 2006 Neil D. Lawrence



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
  set(handle, 'CData', reshape(imageValues, imageSize(1), imageSize(2))');
else
  set(handle, 'CData', reshape(imageValues, imageSize(1), imageSize(2)));
end
