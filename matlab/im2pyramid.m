function Yout = im2pyramid(Yin, height, width, layers)

nSamples = size(Yin,1);
dOld = size(Yin,2);
nH = height;
nW = width;
dNew = 0;
for l=1:layers
    dNew = dNew + nH*nW;
    nH = nH/2;
    nW = nW/2;
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
        idx = idx+newSize+1;
    end    
end