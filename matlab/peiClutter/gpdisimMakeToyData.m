function [genes, params] = gpdisimMakeToyData(numGenes, numPoints)

% GPDISIMMAKETOYDATA Generate a toy data set for the DISIM model.
  
  
t = (1:numPoints)';

delta = .1;
sigma = gamrnd(1, 10);
D = gamrnd(ones(1, numGenes), 10);
S = gamrnd(ones(1, numGenes), 10);
B = gamrnd(ones(1, numGenes), 10);
l = gamrnd(1, 1);

kernType1 = {'multi'};
tieWidth = [1]; % These are the indices of the inverse widths which
                % need to be constrained to be equal.
kernType1{2} = 'rbf';
for i = 1:numGenes
  kernType1{i+2} = 'disim';
  if i==1
    tieDelta = [3];
    tieWidth = [tieWidth, 4];
    tieSigma = [5];
  end
  if i>1
    tieDelta = [tieDelta tieDelta(end)+5];
    tieWidth = [tieWidth tieWidth(end)+5];
    tieSigma = [tieSigma tieSigma(end)+5];
  end
end
tieParam = {tieDelta, tieWidth, tieSigma};


k = kernCreate(t, kernType1);
k = modelTieParam(k, tieParam);
params = kernExtractParam(k);

params(1) = log(l);
params(2) = 0;
params(3) = log(delta);
params(4) = log(sigma);
params(5:2:end) = log(D);
params(6:2:end) = log(S);

mu = B./D;

k = kernExpandParam(k, params);
c = kernCompute(k, t);
m = [ ];

ind = 1:numPoints;
m(ind) = 2;
ind = ind + numPoints;
for i = 1:numGenes
  m(ind) = mu(i)*ones(numPoints, 1);
  ind = ind + numPoints;
end

x = mvnrnd(m, c)';
genes = x + .1 * randn(size(x)); %[zeros(numPoints, 1); randn(numPoints*numGenes, 1)];
genes = reshape(genes, [numPoints, numGenes+1]);
