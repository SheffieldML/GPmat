function doubleMatrixWriteToFID(val, FID)
  
% DOUBLEMATRIXWRITETOFID Writes a double matrix to an FID.
% FORMAT
% DESC writes a double matrix to a stream.
% ARG val : matrix to place in file.
% ARG FID : stream to write to.
%
% SEEALSO : doubleMatrixReadFromFID, matrixWriteToFID
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% MLTOOLS
  
writeVersionToFID(FID, 0.2);
writeStringToFID(FID, 'baseType', 'matrix');
writeStringToFID(FID, 'type', 'doubleMatrix');
writeIntToFID(FID, 'numRows', size(val, 1));
writeIntToFID(FID, 'numCols', size(val, 2));
for i = 1:size(val, 1)
  for j = 1:size(val, 2)-1
    fprintf(FID, '%1.17e ', val(i, j));
  end
  fprintf(FID, '%1.17e\n', val(i, end));
end
  
