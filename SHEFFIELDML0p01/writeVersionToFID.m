function writeVersionToFID(FID, val)

% WRITEVERSIONTOFID Writes a version to an FID.
%
%	Description:
%
%	WRITEVERSIONTOFID(FID, VAL) writes a version from a stream.
%	 Arguments:
%	  FID - stream to write to.
%	  VAL - value of version to place in file.
%	
%
%	See also
%	READVERSIONFROMFID, WRITESTRINGTOFID


%	Copyright (c) 2008 Neil D. Lawrence

  
writeStringToFID(FID, 'version', num2str(val));
  
