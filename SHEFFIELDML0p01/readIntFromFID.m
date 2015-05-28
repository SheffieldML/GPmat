function val = readIntFromFID(FID, string)

% READINTFROMFID Read an integer from an FID.
%
%	Description:
%
%	VAL = READINTFROMFID(FID, NAME) reads an integer from a stream.
%	 Returns:
%	  VAL - value of variable in file.
%	 Arguments:
%	  FID - stream to read from.
%	  NAME - name of integer.
%	readBoolFromFID, readStringFromFID
%	
%
%	See also
%	WRITEINTTOFID, READDOUBLEFROMFID, READINTFROMFID, 


%	Copyright (c) 2008 Neil D. Lawrence

  
  
val = str2num(readStringFromFID(FID, string));
