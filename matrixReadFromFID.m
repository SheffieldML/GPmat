function X = matrixReadFromFID(FID, varargin)
  
% MATRIXREADFROMFID Read a matrix from an FID.
% FORMAT
% DESC reads a matrix from an FID.
% ARG FID: the file ID to read the matrix from.
% RETURN X : the returned matrix read from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% SEEALSO : modelReadFromFID, doublematrixReadFromFID

% MLTOOLS

modelType = readStringFromFID(FID, 'type');
feval = str2func([modelType 'ReadFromFID']);
X = feval(FID, varargin{:});
