function lineStr = getline(FID, comment)

% GETLINE Get a line from a file, but ignore it if it starts with a comment (default #).

% NDLUTIL


if nargin < 2
  comment = '#';
end
if length(comment)~=1
  error('Comment can only be one character in length');
end
lineStr = fgetl(FID);
if length(lineStr)==0
  return
end
while lineStr(1)==comment
  lineStr = fgetl(FID);
end