function noise = noiseTest(noiseType);

% NOISETEST Run some tests on the specified noise model.

if iscell(noiseType)
  % compound noise type
  noise.type = 'cmpnd';
  for i = 1:length(noiseType)
    noise.comp{i}.type = noiseType{i};
  end
else
  noise.type = noiseType;
end
noise.C = 10;
noise.numProcess = 3;
numData = 10;
noise = noiseParamInit(noise);
%noise.eta = 1e-10;% 1/(2*noise.C);
mu = randn(numData, noise.numProcess);
varsigma = randn(numData, noise.numProcess).^2;
y = noiseOut(noise, mu, varsigma);
epsilon = 1e-6;
fprintf('y values\n');
disp(y)
params = noiseExtractParam(noise);
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  noise = noiseExpandParam(params, noise);
  Lplus(i) = noiseLogLikelihood(noise, mu, varsigma, y);
  params(i) = origParams(i) - epsilon;
  noise = noiseExpandParam(params, noise);
  Lminus(i) = noiseLogLikelihood(noise, mu, varsigma, y);
end
params = origParams;
noise = noiseExpandParam(params, noise);
[void, names] = noiseExtractParam(noise);
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = noiseGradientParam(noise, mu, varsigma, y);

paramMaxDiff = max(max(abs(gLDiff-g)));
% l = 0;
% for i = 1:length(names)
%   if l < length(names{i})
%     l = length(names{i});
%   end
% end

% fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
% for i = 1:length(names)
%   spaceLen = l - length(names{i});
%   space = char(repmat(32, 1, spaceLen));
%   fprintf([space names{i} ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
%            g(i), gLDiff(i), gLDiff(i) - g(i));
% end
Lplus = zeros(size(mu));
Lminus = zeros(size(mu));
origMu = mu;
for i = 1:size(mu, 1)
  for j = 1:size(mu, 2)
    mu = origMu;
    mu(i, j) = origMu(i, j) + epsilon;
    Lplus(i, j) = noiseLogLikelihood(noise, mu, varsigma, y);
    mu(i, j) = origMu(i, j) - epsilon;
    Lminus(i, j) = noiseLogLikelihood(noise, mu, varsigma, y);
  end
end
mu = origMu;
gMuDiff = .5*(Lplus - Lminus)/epsilon;

Lplus = zeros(size(varsigma));
Lminus = zeros(size(varsigma));
origVarsigma = varsigma;
for i = 1:size(varsigma, 1)
  for j = 1:size(varsigma, 2)
    varsigma = origVarsigma;
    varsigma(i, j) = origVarsigma(i, j) + epsilon;
    Lplus(i, j) = noiseLogLikelihood(noise, mu, varsigma, y);
    varsigma(i, j) = origVarsigma(i, j) - epsilon;
    Lminus(i, j) = noiseLogLikelihood(noise, mu, varsigma, y);
  end
end
varsigma = origVarsigma;
gVarsigmaDiff = .5*(Lplus - Lminus)/epsilon;



[g, gvs] = noiseGradVals(noise, mu, varsigma, y);

vsMaxDiff = max(max(abs(gvs-gVarsigmaDiff)));
muMaxDiff = max(max(abs(g-gMuDiff)));


fprintf('Param max diff: %2.6f.\n', paramMaxDiff)
fprintf('Mu max diff: %2.6f.\n', muMaxDiff)
fprintf('Varsigma max diff: %2.6f.\n', vsMaxDiff)