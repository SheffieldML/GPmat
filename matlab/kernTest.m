function kernRet = kernTest(kernType, numIn, tieParamNames, flags);

% KERNTEST Run some tests on the specified kernel.
% FORMAT
% DESC runs some tests on a kernel with the specified type to ensure it is
% correctly implemented.
% ARG kernType : type of kernel to test. For example, 'rbf' or
% {'cmpnd', 'rbf', 'lin', 'white'}.
% ARG numIn : the number of input dimensions (default is 4).
% ARG tieParamNames : cell array of regular expressions for parameter
% names that should be tied (default is none).
% ARG flags : vector with a series of binary flags indicating the type of
% the kernel used. Currently there are only two flags indicating: (1) the
% normalisation of the kernel; (2) its stationarity. Note that for some
% kernels these flags may be pointless (e.g. there is not a normalised
% version of the white kernel and it is always stationary). In these cases
% the kernel computation should simply ignore the flags. By default they
% are set to false.
% RETURN kern : the kernel that was generated for the tests.
% 
% FORMAT
% DESC runs some tests on a given kernel structure to ensure it is
% correctly implemented.
% ARG kern : kernel structure to test.
% ARG numIn : the number of input dimensions (default is 4).
% ARG tieParamNames : cell array of regular expressions for parameter
% names that should be tied (default is none).
% ARG flags : vector with a series of binary flags indicating the type of
% the kernel used. Currently there are only two flags indicating: (1) the
% normalisation of the kernel; (2) its stationarity. Note that for some
% kernels these flags may be pointless (e.g. there is not a normalised
% version of the white kernel and it is always stationary). In these cases
% the kernel computation should simply ignore the flags. By default they
% are set to false.
% RETURN kern : the kernel as it was used in the tests.
% 
% SEEALSO : kernCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2007
%
% MODIFICATIONS : Antti Honkela, 2007
%
% MODIFICATIONS : David Luengo, 2009

% KERN

if nargin < 2
  numIn = 4;
end
if nargin < 3
  tieParamNames = {};
end
if nargin < 4
    isNormalised = false;
    isStationary = false;
else
    isNormalised = flags(1);
    isStationary = flags(2);
end
numData = 20;

if isstruct(kernType)
  kern = kernType;
  kernType = kern.type;
  [params, paramnames] = kernExtractParam(kern);
  paramExpand = eye(length(params));
  paramPack = paramExpand;
  toRemove = [];
else
  kern = kernCreate(numIn, kernType);
  if exist([kern.type 'KernSetIndex'])==2 
    for i = 1:length(kern.comp)
      if rand(1)>0.5
        indices = randperm(numIn);
        indices = indices(1:ceil(rand(1)*numIn));
        kern = kernSetIndex(kern, i, indices);
      end
    end
  end
  if isfield(kern, 'positiveTime') && kern.positiveTime 
    % For convolutional kernels starting at t=0 it does not make sense to use
    % negative inputs...
    x = abs(randn(numData, numIn));
    x2 = abs(randn(numData/2, numIn));
  elseif isfield(kern, 'indexValues') && kern.indexValues
    % For index covariance function.
    x = kern.indices(randi(length(kern.indices), numData, 1))';
    x2 = kern.indices(randi(length(kern.indices), numData, 1))';
  else
    x = randn(numData, numIn);
    x2 = randn(numData/2, numIn);
  end
  

  % Set the parameters randomly.
  [params, paramnames] = kernExtractParam(kern);
  if iscell(tieParamNames),
    tieParams = cell(size(tieParamNames));
    paramExpand = eye(length(params));
    toRemove = [];
    for l = 1:length(tieParamNames),
      %ties = strfind(paramnames, tieParamNames{l});
      ties = regexp(paramnames, tieParamNames{l});
      tieParams{l} = [];
      for k = 1:length(ties),
	if ~isempty(ties{k}),
	  tieParams{l} = [tieParams{l}, k];
	end
      end
      if ~isempty(tieParams{l}),
	paramExpand(:, tieParams{l}(1)) = sum(paramExpand(:, tieParams{l}), 2);
	toRemove = [toRemove, tieParams{l}(2:end)];
      end
    end
    paramExpand(:, sort(toRemove)) = [];
  else
    paramExpand = tieParamNames;
    toRemove = [];
    for k=1:size(paramExpand, 2),
      I = find(paramExpand(:, k));
      toRemove = [toRemove, I(2:end)'];
    end
  end
  paramPack = paramExpand' ./ repmat(sum(paramExpand', 2), [1, size(paramExpand, 1)]);
  params = params * paramPack';
  params = randn(size(params))./sqrt(randn(size(params)).^2);
  kern = kernExpandParam(kern, params * paramExpand');
end
% Setting the normalisation and stationarity of the kernel according to the
% flags.
if isfield(kern, 'isNormalised')
    kern.isNormalised = isNormalised;
end
if isfield(kern, 'isStationary')
    kern.isStationary = isStationary;
end

% Test for positive definiteness
K = kernCompute(kern, x);
e = eig(K);
if min(e) > 0,
  fprintf('The kernel is positive definite.\n');
else
  fprintf('The kernel is not positive definite: max eig %g, min eig %g\n', ...
	  max(e), min(e));
end

covGrad = ones(size(K));
epsilon = 1e-6;
params = kernExtractParam(kern) * paramPack';
origParams = params;
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params * paramExpand');
  Lplus(i) = full(sum(sum(kernCompute(kern, x))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params * paramExpand');
  Lminus(i) = full(sum(sum(kernCompute(kern, x))));
end
params = origParams;
kern = kernExpandParam(kern, params * paramExpand');
[void, names] = kernExtractParam(kern);
names(toRemove) = [];
gLDiff = .5*(Lplus - Lminus)/epsilon;
g = kernGradient(kern, x, covGrad) * paramExpand;


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
    fprintf([space names{i} ':\t%4.6g\t%4.6g\t%4.6g\n'], ...
            g(i), gLDiff(i), gLDiff(i) - g(i));
  end
end
try 
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
catch
  fprintf('kernGradX has an error.\n')
  warning(lasterr)
  xMaxDiff = 0;
  xDiagMaxDiff = 0;
end

K = kernCompute(kern, x);
traceK =  full(trace(K));
K2 = kernDiagCompute(kern, x);
traceK2 = full(sum(K2));
traceDiff = traceK - traceK2; 
%if abs(traceDiff) > 2*epsilon,
%  fprintf('kernDiagCompute is not in sync with kernCompute.\n')
%  fprintf('diag(kernCompute)\tkernDiagCompute')
%  disp([diag(K), K2])
%end

covGrad = ones(size(kernCompute(kern, x, x2)));
epsilon = 1e-6;
params = kernExtractParam(kern) * paramPack';
origParams = params;
Lplus = zeros(size(params));
Lminus = zeros(size(params));
for i = 1:length(params);
  params = origParams;
  params(i) = origParams(i) + epsilon;
  kern = kernExpandParam(kern, params * paramExpand');
  Lplus(i) = full(sum(sum(kernCompute(kern, x, x2))));
  params(i) = origParams(i) - epsilon;
  kern = kernExpandParam(kern, params * paramExpand');
  Lminus(i) = full(sum(sum(kernCompute(kern, x, x2))));
end
params = origParams;
kern = kernExpandParam(kern, params * paramExpand');
[void, names] = kernExtractParam(kern);
names(toRemove) = [];
gL2Diff = .5*(Lplus - Lminus)/epsilon;
g = kernGradient(kern, x, x2, covGrad) * paramExpand;

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
    fprintf([space names{i} ':\t%4.6g\t%4.6g\t%4.6g\n'], ...
            g(i), gL2Diff(i), gL2Diff(i) - g(i));
  end
  pause(0);
end


fprintf('Trace max diff: %2.6g.\n', traceDiff);
fprintf('Param max diff: %2.6g.\n', paramMaxDiff)
fprintf('Param X2 max diff: %2.6g.\n', param2MaxDiff)
fprintf('X max diff: %2.6g.\n', xMaxDiff)
fprintf('XDiag max diff: %2.6g.\n', xDiagMaxDiff)
fprintf('\n');

if nargout > 0
  kernRet = kern;
else
  kernDisplay(kern);
end

% We don't test kernCompute(kern, x, x2) here at all!
