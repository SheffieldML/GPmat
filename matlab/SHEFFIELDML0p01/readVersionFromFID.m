function val = readVersionFromFID(FID)

% READVERSIONFROMFID Read version number from an FID.
%
%	Description:
%
%	VAL = READVERSIONFROMFID(FID) reads version number from a stream.
%	 Returns:
%	  VAL - value of variable in file.
%	 Arguments:
%	  FID - stream to read from.
%	
%
%	See also
%	WRITEVERSIONTOFID, READDOUBLEFROMFID, READSTRINGFROMFID


%	Copyright (c) 2008 Neil D. Lawrence

  
val = str2num(readStringFromFID(FID, 'version'));
