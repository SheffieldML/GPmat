function writeStringToFID(FID, name, val)

% WRITESTRINGTOFID Writes a string to an FID.
%
%	Description:
%
%	WRITESTRINGTOFID(FID, NAME, VAL) writes an string from a stream.
%	 Arguments:
%	  FID - stream to write from.
%	  NAME - name of string.
%	  VAL - value of variable to place in file.
%	
%
%	See also
%	READSTRINGFROMFID, WRITEINTTOFID, WRITEBOOLTOFID, WRITEDOUBLETOFID


%	Copyright (c) 2008 Neil D. Lawrence

  
fprintf(FID, [name '=' val '\n']);
