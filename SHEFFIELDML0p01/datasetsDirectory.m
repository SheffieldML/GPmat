function dir = datasetsDirectory

% DATASETSDIRECTORY Returns directory where data is stored.
%
%	Description:
%
%	DIR = DATASETSDIRECTORY returns the directory where this file is
%	located, it is assumed that any data sets downloaded have this
%	directory as their base directory.
%	 Returns:
%	  DIR - the directory where this m-file is stored.
%	
%
%	See also
%	LVMLOADDATA, MAPLOADDATA


%	Copyright (c) 2005, 2006 Neil D. Lawrence



% by default return the directory where this file is.
fullSpec = which('datasetsDirectory');
ind = max(find(fullSpec == filesep));
dir = fullSpec(1:ind);
