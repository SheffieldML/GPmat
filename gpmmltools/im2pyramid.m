function Yout = im2pyramid(Yin, height, width, layers)
% Yout = im2pyramid(Yin, height, width, layers)
%            Expands a data matrix that has been generated from images by 
%            concatenating it with the pyramid representation of each
%            image (using a low-pass filter)
% FORMAT
% Yin: Input design matrix, each row is a sample (built by vectorysing an
%      image).
% height: Hight of the original images.
% width: Width of the original images.
% layers: Number of layers used to built the image pyramids.
%
% Note: requires MatLab's Image Processing toolbox.
%
% COPYRIGHT: Teo de Campos 2013
%

% SHEFFIELDML

nSamples = size(Yin,1);
dOld = size(Yin,2);
nH = height;
nW = width;
dNew = 0;
for l=1:layers
    dNew = dNew + nH*nW;
    nH = ceil(nH/2);
    nW = ceil(nW/2);
end

Yout = zeros(nSamples, dNew);
Yout(:,1:size(Yin,2)) = Yin;
for s = 1:nSamples,
    idx = dOld;
    img = reshape(Yin(s,:), height, width);
    for l=2:layers
        img = impyramid(img, 'reduce');
        newSize = numel(img);
        Yout(s,idx+1:idx+newSize) = img(:)';
        idx = idx+newSize;
    end    
end