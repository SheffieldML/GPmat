function model = modelLoadResult(type, dataSet, number)

% MODELLOADRESULT Load a previously saved result.
% FORMAT
% DESC loads a previously saved model result.
% ARG type : the type of model to load.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the model run to load.
% RETURN model : the saved model.
%
% SEEALSO : modelLoadResult
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% MLTOOLS
  
  origDataSet = dataSet;
  dataSet(1) = upper(dataSet(1));
  origType = type;
  type(1) = upper(type(1));
  fileName = ['dem' dataSet type num2str(number)];

  fhandle = [type 'LoadResult'];
  if exist(fhandle)==2
    % There is load result code, use it.
    fhandle = str2func(fhandle);
    model = fhandle(origDataSet, number);
  else
    fhandle = [origType 'Reconstruct'];
    if exist(fhandle)==2
      % There is reconstruct code, use it to reconstruct.
      [Y, lbls] = lvmLoadData(origDataSet);
      load(fileName);
      varName = [origType 'Info'];
      eval(['model = ' fhandle '(' varName ', Y);'])
    else
      % No code is provided, just load the file.
      load(fileName);
    end
  end
end