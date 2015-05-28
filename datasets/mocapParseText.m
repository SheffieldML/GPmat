function [points, pointNames, times] = mocapParseText(fileName);

% MOCAPPARSETEXT Parse a motion capture text file.
% FORMAT
% DESC parses a motion capture text file and returns point
% positions, point names and times.
% ARG fileName : the file name in which the data is stored.
% RETURN points : point positions extracted from the file.
% RETURN pointNames : the names of the points in the files.
% RETURN times : the timing of each individual frame in the file.
%
% SEEALSO : mocapConnections, mocapLoadTextData
%
% COPYRIGHT : Neil D. Lawrence, 2006

% DATASETS


fid = fopen(fileName);
readLine = fgets(fid);

[token, rem] = strtok(readLine);
[token, rem] = strtok(rem); % skipping the 'field' field
[token, rem] = strtok(rem); % skipping the 'time' field

i = 1;
while (~isempty(token))
    pointNames{i}= token(1:end-2);
    i = i+1;
    [token, rem] = strtok(rem);
    [token, rem] = strtok(rem);
    [token, rem] = strtok(rem);
end 

numPoints = length(pointNames);
S = fscanf(fid, '%f ');
numData = length(S)/(numPoints*3+2);
S = reshape(S, numPoints*3+2, numData)';
field = S(:, 1);
times = S(:, 2);
S(find(S==-9999.99)) = NaN;
points{1} = S(:, 3:3:end);
points{2} = S(:, 4:3:end);
points{3} = S(:, 5:3:end);

fclose(fid);
