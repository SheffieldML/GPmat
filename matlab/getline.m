function lineStr = getline(FID)

% GETLINE Get a line from a file, but ignore it if it starts with #.

% NDLUTIL

lineStr = fgetl(FID);
if length(lineStr)==0
  return
end
while lineStr(1)=='#' 
  lineStr = fgetl(FID);
end