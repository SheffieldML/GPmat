function kernRet = kernTest(kernType, numIn);

% KERNTEST Run some tests on the specified kernel.
% FORMAT
% DESC runs some tests on a kernel with the specified type to ensure it is
% correctly implemented.
% ARG kernType : type of kernel to test. For example, 'rbf' or
% {'cmpnd', 'rbf', 'lin', 'white'}.
% RETURN kern : the kernel that was generated for the tests.
% 
% FORMAT
% DESC runs some tests on a given kernel structure to ensure it is
% correctly implemented.
% ARG kern : kernel structure to test.
% ARG numIn : the number of input dimensions.
% RETURN kern : the kernel as it was used in the tests.
% 
% SEEALSO : kernCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2007

% KERN

if nargin < 2
  numIn = 4;
end
numData = 20;

% Generate some x positions.
x = randn(numData, numIn);
x2 = randn(numData/2, numIn);
if isstruct(kernType)
  kern = kernType;
else
  kern = kernCreate(x, kernType);
  if exist([kern.type 'KernSetIndex'])==2 
    for i = 1:length(kern.comp)
      if rand(1)>0.5
        indices = randperm(numIn);
        indices = indices(1:ceil(rand(1)*numIn));
        kern = kernSetIndex(kern, i, indices);
      end
    end
  end
  
  % Set the parameters randomly.
  params = kernExtractParam(kern);
  params = randn(size(params))./sqrt(randn(size(params)).^2);
  kern = kernExpandParam(kern, params);
end

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

covGrad = ones(numData, numData/2);
epsilon = 1e-6;
params = kernExtractParam(kern);
origParams = params;
Lplus = zeros(size(params));
Lminus = zeros(size(params));
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params);
  Lplus(i) = full(sum(sum(kernCompute(kern, x, x2))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params);
  Lminus(i) = full(sum(sum(kernCompute(kern, x, x2))));
end
params = origParams;
kern = kernExpandParam(kern, params);
[void, names] = kernExtractParam(kern);
gL2Diff = .5*(Lplus - Lminus)/epsilon;
g = kernGradient(kern, x, x2, covGrad);

param2MaxDiff = max(max(abs(gL2Diff-g)));
if param2MaxDiff > 2*epsilon
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
            g(i), gL2Diff(i), gL2Diff(i) - g(i));
  end
end


fprintf('Trace max diff: %2.6f.\n', traceDiff);
fprintf('Param max diff: %2.6f.\n', paramMaxDiff)
fprintf('Param X2 max diff: %2.6f.\n', param2MaxDiff)
fprintf('X max diff: %2.6f.\n', xMaxDiff)
fprintf('XDiag max diff: %2.6f.\n', xDiagMaxDiff)
fprintf('\n');

if nargout > 0
  kernRet = kern;
else
  kernDisplay(kern);
end

% We don't test kernCompute(kern, x, x2) here at all!