function [model, origI] = ivmReadFromFID(FID)

% IVMREADFROMFID Load from a FID produced by the C++ implementation.
% FORMAT
% DESC loads in from a file stream the data format produced by the
% C++ implementation of the IVM.
% ARG FID : the file ID from where the data is loaded.
% RETURN model : the model loaded in from the file.
% model.
% RETURN I : the indices of the active points in the original data set.
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% SEEALSO : ivmReadFromFile

% IVM

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'ivmVersion'))
  error('Incorrect file format.')
end
if(~strcmp(tokens{2}, '0.1'))
  error('Incorrect file version.')
end
version = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'activeSetSize'))
  error('Incorrect file format.')
end
activeSetSize = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numProcesses'))
  error('Incorrect file format.')
end
numProcesses = str2num(tokens{2});

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'numFeatures'))
  error('Incorrect file format.')
end
numFeatures = str2num(tokens{2});
kern = kernReadFromFID(FID);
noise = noiseReadFromFID(FID);

X = zeros(activeSetSize, numFeatures);
y = zeros(activeSetSize, numProcesses);
m = zeros(activeSetSize, numProcesses);
beta = zeros(activeSetSize, numProcesses);

origI = zeros(activeSetSize, 1);

for i=1:activeSetSize
  lineStr = getline(FID);
  tokens = tokenise(lineStr);
  for j = 1
    origI(i) = str2num(tokens{j});
  end
  for j=1:numProcesses
    y(i, j) = str2num(tokens{1+j});
  end
  for j=1:numProcesses
    m(i, j) = str2num(tokens{1+numProcesses+j});
  end
  for j=1:numProcesses
    beta(i, j) = str2num(tokens{1+2*numProcesses+j});
  end
  for j=numProcesses*3+2:length(tokens)
    str = tokens{j};
    ind = find(str==':');
    % TODO -- check that ':' is in the string.
    featStr = str(1:ind-1);
    featValStr = str(ind+1:end);
    featNum = str2num(featStr);
    featVal = str2num(featValStr);
    X(i, featNum) = featVal;
  end
end

ivmInfo.J = [];
ivmInfo.I = 1:activeSetSize;
ivmInfo.m = m;
ivmInfo.beta = beta;

model = ivmReconstruct(kern, noise, ivmInfo, X, y);
