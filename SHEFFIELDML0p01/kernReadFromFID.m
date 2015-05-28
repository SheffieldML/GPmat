function kern = kernReadFromFID(FID, version)

% KERNREADFROMFID Load from an FID written by the C++ implementation.
%
%	Description:
%
%	KERN = KERNREADFROMFID(FID) loads in from a file stream the data
%	format produced by C++ implementations.
%	 Returns:
%	  KERN - the kernel loaded in from the file.
%	 Arguments:
%	  FID - the file ID from where the data is loaded.
%	
%
%	See also
%	MODELREADFROMFID, KERNCREATE, KERNREADPARAMSFROMFID


%	Copyright (c) 2005, 2006, 2008, 2009 Neil D. Lawrence


  if nargin < 2
    version = 0.2;
  end
  if version <= 0.11
    version = readDoubleFromFID(FID, 'kernVersion');
  end
  type = readStringFromFID(FID, 'type');
  kern = kernCreate(zeros(1), type);

  kern = kernReadParamsFromFID(kern, FID, version);
