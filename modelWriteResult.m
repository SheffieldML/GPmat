function fileName = modelWriteResult(model, dataset, number)

% MODELWRITERESULT Write a model to file.
% FORMAT
% DESC writes a model to file.
% ARG model : the model to be written to file.
% ARG datasetName : data set name.
% ARG number : the experiment number.
% RETURN fileName : the file name used.
%
% SEEALSO : modelCreate, modelLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS
  origDataset = dataset;
  dataset(1) = upper(dataset(1));
  type = model.type;
  type(1) = upper(type(1));
  fileName = ['dem' dataset type num2str(number)];
  if isoctave
    fileName = [fileName '.mat'];
  end
  fhandle = [model.type 'WriteResult'];
  if exist(fhandle)==2
    % There is write file code, use it.
    fhandle = str2func(fhandle);
    fileName = fhandle(model, origDataset, number);
  else
    fhandle = [model.type 'Deconstruct'];
    if exist(fhandle)==2
      % There is deconstruct code, use it to deconstruct.
      varName = [model.type 'Info'];
      eval([varName ' = ' fhandle '(model);'])
      save(fileName, varName);
    else
      % No code is provided, just save model.
      save(fileName, 'model');
    end
  end
end
