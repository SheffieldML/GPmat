function kern = kernTest(kernType);

% KERNTEST Run some tests on the specified kernel.

% KERN


numData = 20;
numIn = 2;

% Generate some x positions.
x = randn(numData, numIn);
kern = kernel(x, kernType);
kern = kernParamInit(kern);

% Set the parameters randomly.
params = kernExtractParam(kern);
params = randn(size(params))./sqrt(randn(size(params)).^2);
kern = kernExpandParam(kern, params);

covGrad = ones(numData);
epsilon = 1e-6;
params = kernExtractParam(kern);
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params);
  Lplus(i) = full(sum(sum(kernCompute(kern, x))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params);
  Lminus(i) = full(sum(sum(kernCompute(kern, x))));
end
params = origParams;
kern = kernExpandParam(kern, params);
[void, names] = kernExtractParam(kern);
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = kernGradient(kern, x, covGrad);

paramMaxDiff = max(max(abs(gLDiff-g)));
if paramMaxDiff > 2*epsilon
  l = 0;
  for i = 1:length(names)
    if l < length(names{i})
      l = length(names{i});
    end
  end
  
  fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
  for i = 1:length(names)
    spaceLen = l - length(names{i});
    space = char(repmat(32, 1, spaceLen));
    fprintf([space names{i} ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
            g(i), gLDiff(i), gLDiff(i) - g(i));
  end
end

Lplus = zeros(size(x));
Lminus = zeros(size(x));
gx = zeros(size(x));
origX = x;
for i = 1:size(x, 1)
  for j = 1:size(x, 2)
    x = origX;
    x(i, j) = origX(i, j) + epsilon;
    K = kernCompute(kern, x);
    Lplus(i, j) =  full(sum(sum(K)));
    LplusDiag(i, j) = full(trace(K));
    x(i, j) = origX(i, j) - epsilon;
    K = kernCompute(kern, x);
    Lminus(i, j) = full(sum(sum(K)));
    LminusDiag(i, j) = full(trace(K));
  end
  x = origX;
  gx(i, :) = 2*sum(kernGradX(kern, x(i, :), x), 1);
  gxDiag(i, :) = kernDiagGradX(kern, x(i, :));
end

gXDiff = .5*(Lplus - Lminus)/epsilon;
xMaxDiff = max(max(abs(gx-gXDiff)));

if xMaxDiff > 2*epsilon
  fprintf('gX\n')
  disp(gx)
  fprintf('gXDiff\n')
  disp(gXDiff)
end

gXDiagDiff = .5*(LplusDiag - LminusDiag)/epsilon;
xDiagMaxDiff = max(max(abs(gxDiag-gXDiagDiff)));

if xDiagMaxDiff > 2*epsilon
 fprintf('gxDiag\n')
 disp(gxDiag)
 fprintf('gXDiagDiff\n')
 disp(gXDiagDiff)
end

K = kernCompute(kern, x);
traceK =  full(trace(K));
traceK2 = full(sum(kernDiagCompute(kern, x)));
traceDiff = traceK - traceK2; 
  
fprintf('Trace max diff: %2.6f.\n', traceDiff);
fprintf('Param max diff: %2.6f.\n', paramMaxDiff)
fprintf('X max diff: %2.6f.\n', xMaxDiff)
fprintf('XDiag max diff: %2.6f.\n', xDiagMaxDiff)
