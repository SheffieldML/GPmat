function model = gpReadFromFID(FID, varargin)

% GPREADFROMFID Load from a FID produced by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by the
% C++ GP implementation.
% ARG FID : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2008, 2009
%
% SEEALSO : gpReadFromFile

% GP

try
  version = readVersionFromFID(FID);
catch
  if strcmp(lasterr, 'Incorrect file format')
    version = readDoubleFromFID(FID, 'gplvmVersion');
  else
    error(lasterr)
  end
end

if version > 0.2 
  error('Incorrect file version.')
end

if version>0.11
  baseType = readStringFromFID(FID, 'baseType')
    if ~strcmp(baseType, 'dataModel') && ~strcmp(baseType, 'mapmodel')
    error('Incorrect base type in file.')
  end
  type = readStringFromFID(FID, 'type')
  if ~strcmp(type, 'gp')
    error('Incorrect type in file.')
  end
end
numData = readIntFromFID(FID, 'numData');
dataDim = readIntFromFID(FID, 'outputDim');
inputDim = readIntFromFID(FID, 'inputDim');
sparseApprox = readIntFromFID(FID, 'sparseApproximation');
numActive = readIntFromFID(FID, 'numActive');
if sparseApprox
  beta = modelReadFromFID(FID);
  beta = beta(1, 1);
end

learnScale = readBoolFromFID(FID, 'learnScale');
learnBias = readBoolFromFID(FID, 'learnBias');
scale = modelReadFromFID(FID);
bias = modelReadFromFID(FID);

kern = modelReadFromFID(FID);
noise = modelReadFromFID(FID);

if sparseApprox
  fixInducing = readBoolFromFID(FID, 'fixInducing');
  X_u = modelReadFromFID(FID);
end

X = varargin{1};
y = varargin{2};
% X = zeros(numData, inputDim);
% y = zeros(numData, dataDim);

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
 case 4
  approxType = 'dtcvar';
end
options = gpOptions(approxType);
options.numActive = numActive;
options.kern = kern;
if sparseApprox
  options.fixInducing = fixInducing;
end
model = gpCreate(inputDim, size(y, 2), X, y, options);
model.X = X;
if sparseApprox
  model.X_u = X_u;
  model.beta = beta;
end
model.scale = scale;
model.bias = bias;
model.m = gpComputeM(model);

% This forces kernel computation.
initParams = gpExtractParam(model);
model = gpExpandParam(model, initParams);
