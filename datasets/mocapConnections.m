function connect = mocapConnections(fileName, pointNames);

% MOCAPCONNECTIONS Give a connection matrix for the motion capture data.
% FORMAT
% DESC load a connection matrix for txt file based mocap data.
% ARG fileName : the file from which to load the connectivity data.
% ARG pointNames : the names of the points in the file, with the
% ordering matching orderings in any data files to be loaded in (as
% returned by MOCAPPARSETEXT).
%
% SEEALSO : mocapLoadTextData, mocapParseText
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% DATASETS


fid = fopen(fileName);
i = 1;
rem = fgets(fid);	
while(rem ~= -1)		
  [token, rem] = strtok(rem, ',');
  connections{i, 1} = fliplr(deblank(fliplr(deblank(token))));
  [token, rem] = strtok(rem, ',');
  connections{i, 2} = fliplr(deblank(fliplr(deblank(token))));
  i = i + 1;
  rem = fgets(fid);	
end

connect = zeros(length(pointNames));
fclose(fid);
for i = 1:size(connections, 1);
  for j = 1:length(pointNames)
    if strcmp(pointNames{j}, connections{i, 1}) | ...
          strcmp(pointNames{j}, connections{i, 2})
      for k = 1:length(pointNames)
        if k == j
          break
        end
        if strcmp(pointNames{k}, connections{i, 1}) | ...
              strcmp(pointNames{k}, connections{i, 2})
          connect(j, k) = 1;
        end
      end
    end
  end
end
connect = sparse(connect);      
