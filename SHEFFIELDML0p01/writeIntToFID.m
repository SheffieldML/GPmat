function writeIntToFID(FID, name, val)

% WRITEINTTOFID Writes an integer to an FID.
%
%	Description:
%
%	WRITEINTTOFID(FID, NAME, VAL) writes an integer to a stream.
%	 Arguments:
%	  FID - stream to write to.
%	  NAME - name of int.
%	  VAL - value of variable to place in file.
%	
%
%	See also
%	READINTFROMFID, WRITESTRINGTOFID, WRITEBOOLTOFID, WRITEDOUBLETOFID


%	Copyright (c) 2008 Neil D. Lawrence

  
writeStringToFID(FID, name, num2str(val));
