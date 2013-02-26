function [model, origI] = ivmReadFromFID(FID, varargin)

% IVMREADFROMFID Load from a FID produced by the C++ implementation.
%
%	Description:
%
%	[MODEL, I] = IVMREADFROMFID(FID) loads in from a file stream the
%	data format produced by the C++ implementation of the IVM.
%	 Returns:
%	  MODEL - the model loaded in from the file. model.
%	  I - the indices of the active points in the original data set.
%	 Arguments:
%	  FID - the file ID from where the data is loaded.
%	
%
%	See also
%	IVMREADFROMFILE


%	Copyright (c) 2007, 2008 Neil D. Lawrence


numData = readIntFromFID(FID, 'numData');
numProcesses = readIntFromFID(FID, 'outputDim');
numFeatures = readIntFromFID(FID, 'inputDim');
activeSetSize = readIntFromFID(FID, 'activeSetSize');

kern = modelReadFromFID(FID);
noise = modelReadFromFID(FID);

X = zeros(activeSetSize, numFeatures);
y = zeros(activeSetSize, numProcesses);
m = zeros(activeSetSize, numProcesses);
beta = zeros(activeSetSize, numProcesses);

origI = zeros(activeSetSize, 1);

lineStr = getline(FID);
tokens = tokenise(lineStr, '=');
if(length(tokens)~=2 | ~strcmp(tokens{1}, 'activeSet'))
  error('Incorrect file format.')
end
tokens = tokenise(tokens{2}, ' ');
if(length(tokens)~=activeSetSize)
  error('Incorrect file format.')
end
for i = 1:activeSetSize
  origI(i) = str2num(tokens{i});
end
y = modelReadFromFID(FID);
X = modelReadFromFID(FID);
m = modelReadFromFID(FID);
beta = modelReadFromFID(FID);
ivmInfo.J = [];
ivmInfo.I = 1:activeSetSize;
ivmInfo.m = m;
ivmInfo.beta = beta;

model = ivmReconstruct(kern, noise, ivmInfo, X, y);
