function handle = imageVisualise(imageVals, imageSize)

% IMAGEVISUALISE Helper code for showing an image during 2-D visualisation.

handle = imagesc(reshape(imageVals, imageSize(1), imageSize(2))');
