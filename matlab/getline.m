function lineStr = getline(FID, comment)

% GETLINE Get a line from a file.
% FORMAT
% DESC gets a line from a file, but ignore it if it starts with a
% comment character.
% ARG fid : the identity of the file to load from.
% ARG comment : the character that indicates a line is a comment
% (default #).
% RETURN lineStr : the next line in the file.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : fgetl, fopen

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
