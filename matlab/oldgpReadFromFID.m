function [model, lbls] = oldgpReadFromFID(FID)

% GPREADFROMFID Load from a FID produced by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by the
% C++ GP implementation.
% ARG FID : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : gpReadFromFile

% GP

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'gpVersion'))
  error('Incorrect file format.')
end
if(~strcmp(tokens{2}, '0.1'))
  error('Incorrect file version.')
end
version = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'sparseApproximation'))
  error('Incorrect file format.')
end
sparseApprox = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numActive'))
  error('Incorrect file format.')
end
numActive = str2num(tokens{2});

if sparseApprox
  lineStr = getline(FID);
  tokens=tokenise(lineStr, '=');
  if(length(tokens)~=2 | ~strcmp(tokens{1}, 'beta'))
    error('Incorrect file format.')
  end
  beta = str2num(tokens{2});
end

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numData'))
  error('Incorrect file format.')
end
numData = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numProcesses'))
  error('Incorrect file format.')
end
dataDim = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'inputDim'))
  error('Incorrect file format.')
end
inputDim = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, ' ');
for i=1:length(tokens)
  scale(1, i) = str2num(tokens{i});
end

lineStr = getline(FID);
tokens = tokenise(lineStr, ' ');
for i=1:length(tokens)
  bias(1, i) = str2num(tokens{i});
end


kern = kernReadFromFID(FID);
noise = noiseReadFromFID(FID);

if sparseApprox
  lineStr = getline(FID);
  tokens = tokenise(lineStr, ':');
  if strcmp(tokens{1}, 'X_u')
    if str2num(tokens{2})~=inputDim
      error('Incorrect file format.');
    end
  else
    error('Incorrect file format.');
  end
end

X_u = zeros(numActive, inputDim);

for i=1:numActive
  lineStr = getline(FID);
  tokens = tokenise(lineStr);
  for j=1:inputDim
    X_u(i, j) = str2num(tokens{j});
  end
end


lineStr = getline(FID);
tokens = tokenise(lineStr,',');
for i=1:length(tokens)
  subTokens = tokenise(tokens{i}, ':');
  if strcmp(subTokens{1}, 'Y')
    if str2num(subTokens{2})~=dataDim
      error('Incorrect file format.');
    end
  elseif strcmp(subTokens{1}, 'X')
    if str2num(subTokens{2})~=inputDim
      error('Incorrect file format.');
    end
  else
    error('Incorrect file format.');
  end
end

X = zeros(numData, inputDim);
y = zeros(numData, dataDim);

for i=1:numData
  lineStr = getline(FID);
  tokens = tokenise(lineStr);
  for j=1:dataDim
    y(i, j) = str2num(tokens{j});
  end
  for j=1:inputDim
    X(i, j) = str2num(tokens{j+dataDim});
  end
end
warning('Noise model is ignored');

switch sparseApprox
 case 0
  approxType = 'ftc';
 case 1
  approxType = 'dtc';
 case 2
  approxType = 'fitc';
 case 3
  approxType = 'pitc';
end
options = gpOptions(approxType);
options.numActive = numActive;
options.kern = kern;
model = gpCreate(inputDim, size(y, 2), X, y, options);
model.X = X;
model.X_u = X_u;
if sparseApprox
  model.beta = beta;
end
model.scale = scale;
model.bias = bias;
model.m = gpComputeM(model);

% This forces kernel computation.
initParams = gpExtractParam(model);
model = gpExpandParam(model, initParams);
