function lineStr = getline(FID, comment)

% GETLINE Get a line from a file.
%
%	Description:
%
%	LINESTR = GETLINE(FID, COMMENT) gets a line from a file, but ignore
%	it if it starts with a comment character.
%	 Returns:
%	  LINESTR - the next line in the file.
%	 Arguments:
%	  FID - the identity of the file to load from.
%	  COMMENT - the character that indicates a line is a comment
%	   (default #).
%	
%
%	See also
%	FGETL, FOPEN


%	Copyright (c) 2005 Neil D. Lawrence



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
