function imageModify(handle, imageValues, imageSize)

% IMAGEMODIFY Helper code for visualisation of image data.

set(handle, 'CData', reshape(imageValues, imageSize(1), imageSize(2))');
