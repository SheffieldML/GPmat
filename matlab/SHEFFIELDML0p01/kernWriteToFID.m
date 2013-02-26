function kernWriteToFID(kern, FID)

% KERNWRITETOFID Load from an FID written by the C++ implementation.
%
%	Description:
%
%	KERNWRITETOFID(KERN, FID) loads in from a file stream the data
%	format produced by C++ implementations.
%	 Arguments:
%	  KERN - the kernel structure to write to the stream.
%	  FID - the file ID from where the data is loaded.
%	
%
%	See also
%	MODELREADFROMFID, KERNCREATE, KERNREADPARAMSFROMFID


%	Copyright (c) 2005, 2006, 2008 Neil D. Lawrence


writeVersionToFID(FID, 0.2);
writeStringToFID(FID, 'baseType', 'kern');  
writeStringToFID(FID, 'type', kern.type);
kernWriteParamsToFID(kern, FID);
