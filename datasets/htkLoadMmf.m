function [meanValCell, varValCell] = htkLoadMmf(phoneList, fileName)

% HTKLOADMMF File for loading synthesis data from HTK files.
% FORMAT
% DESC loads the means and variances for the phone model from the
% relevant HTK file.
% ARG phoneList : cell array containing the names of the phone you want to load in.
% ARG fileName : the name of the HTK file from which you are loading.
% RETURN meanValCell : the mean values from the file.
% RETURN varValCell : the standard deviation values from the file.
%  
% COPYRIGHT : Neil D. Lawrence, 2009

% DATASETS
  
  inState = false;
  inMean = false;
  inVar = false;
  inHmm = false;
  inStream = false;
  rightPhone = false;
  meanValCell = {};
  varValCell = {};

  fid = fopen(fileName);
  if fid == -1
    error(['No such file name ' fileName])
  end
  readLine = fgets(fid);
  counter = 0;
  data = [];
  while readLine ~= -1
    readLine = fgets(fid);
    counter = counter + 1;
    if strcmp(readLine(1:min(end,10)), '<BEGINHMM>')
      inHmm = true;
      continue
    elseif strcmp(readLine(1:min(end,8)), '<ENDHMM>')
      inHmm = false;
      continue
    elseif strcmp(readLine(1:min(end,2)), '~h')
      splits = stringSplit(readLine, ' ');
      phone = splits{2}(2:end-2);
      match = strcmp(phoneList, phone);
      if any(match)
        i = min(find(match));
        meanValCell{i, 1} = [];
        meanValCell{i, 2} = [];
        meanValCell{i, 3} = [];
        varValCell{i, 1} = [];
        varValCell{i, 2} = [];
        varValCell{i, 3} = [];
        rightPhone = true;
      else
        rightPhone = false;
      end
    end
    if rightPhone &&inHmm 
      if strcmp(readLine(1:min(end,7)), '<STATE>')
        inState = true;
        continue
      end
      if inState && strcmp(readLine(1:min(end,10)), '<STREAM> 1')
        inStream = true;
        continue
      end
      if inStream && strcmp(readLine(1:min(end,6)), '<MEAN>')
        inMean = true;
        continue
      end
      if inStream && strcmp(readLine(1:min(end,10)),'<VARIANCE>')
        inVar = true;
        continue
      end
      if inMean 
        allVals = sscanf(readLine, '%f');
        meanValCell{i, 1} = [meanValCell{i, 1} allVals(1:25)'];
        meanValCell{i, 2} = [meanValCell{i, 2} allVals(26:50)'];
        meanValCell{i, 3} = [meanValCell{i, 3} allVals(51:75)'];
        inMean = false;
        continue
      end
      if inVar
        allVals = sscanf(readLine, '%f');
        varValCell{i, 1} = [varValCell{i, 1} allVals(1:25)'];
        varValCell{i, 2} = [varValCell{i, 2} allVals(26:50)'];
        varValCell{i, 3} = [varValCell{i, 3} allVals(51:75)'];
        inVar = false;
        inState = false;
        inStream = false;
        continue
      end
    end
  end
  fclose(fid);
end
