function X = doubleMatrixReadFromFID(FID)
  
% DOUBLEMATRIXREADFROMFID Read a full matrix from an FID.
% FORMAT
% DESC reads a matrix from an FID.
% ARG FID: the file ID to read the matrix from.
% RETURN X : the returned matrix read from the file.
%
% COPYRIGHT : Neil D. Lawrence, 2008
% 
% SEEALSO : modelReadFromFile

% MLTOOLS
  
numRows = readIntFromFID(FID, 'numRows');
numCols = readIntFromFID(FID, 'numCols');

X = zeros(numRows, numCols);
for i = 1:numRows
  lineStr = getline(FID);
  tokens = tokenise(lineStr, ' ');
  if strcmp(tokens(end),'');
    tokens=tokens(1:end-1);
  end
  if(length(tokens)~=numCols)
    error('Incorrect file format.');
  end
  for j = 1:numCols
    X(i,j) = str2num(tokens{j});
  end
end
