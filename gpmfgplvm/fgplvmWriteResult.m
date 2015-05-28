function fileName = fgplvmWriteResult(model, dataSet, number)

% FGPLVMWRITERESULT Write a FGPLVM result.
% FORMAT
% DESC writess an FGPLVM result.
% ARG model : the model to write.
% ARG dataSet : the name of the data set to write.
% ARG number : the number of the FGPLVM run to write.
% RETURN fileName : the file name used to write.
%
% SEEALSO : fgplvmLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009, 2011
  
% FGPLVM

dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

[kern, noise, fgplvmInfo, X] = fgplvmDeconstruct(model);

save(fileName, 'kern', 'noise', 'fgplvmInfo', 'X');
