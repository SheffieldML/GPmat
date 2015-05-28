function X = doubleMatrixReadFromFID(FID)

% DOUBLEMATRIXREADFROMFID Read a full matrix from an FID.
%
%	Description:
%
%	X = DOUBLEMATRIXREADFROMFID(FID) reads a matrix from an FID.
%	 Returns:
%	  X - the returned matrix read from the file.
%	 Arguments:
%	  FID - the file ID to read the matrix from.
%	
%
%	See also
%	MODELREADFROMFILE


%	Copyright (c) 2008 Neil D. Lawrence

  
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