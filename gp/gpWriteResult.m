function fileName = gpWriteResult(model, dataSet, number)

% GPWRITERESULT Write a GP result.
% FORMAT
% DESC writes a GP result.
% ARG model : the model to write.
% ARG dataSet : the name of the data set to write.
% ARG number : the number of the GP run to write.
% RETURN fileName : the file name used to write.
%
% SEEALSO : gpLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% GP

dataSet(1) = upper(dataSet(1));
type = model.type;
type(1) = upper(type(1));
fileName = ['dem' dataSet type num2str(number)];

[kern, noise, gpInfo] = gpDeconstruct(model);

save(fileName, 'kern', 'noise', 'gpInfo');
