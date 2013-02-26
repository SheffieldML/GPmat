function acclaimPlayFile(fileNameAsf, fileNameAmc, frameLength)

% ACCLAIMPLAYFILE Play motion capture data from a asf and amc file.
%
%	Description:
%
%	ACCLAIMPLAYFILE(FILENAMEASF, FILENAMEAMC, FRAMELENGTH) plays motion
%	capture data by reading in a asf and amc files from disk.
%	 Arguments:
%	  FILENAMEASF - the name of the ASF file to read in.
%	  FILENAMEAMC - the name of the AMC file to read in.
%	  FRAMELENGTH - the length of the frames.
%	
%
%	See also
%	ACCLAIMPLAYFILE, BVHREADFILE, SKELPLAYDATA


%	Copyright (c) 2006 Neil D. Lawrence


skel = acclaimReadSkel(fileNameAsf);
[channels, skel] = acclaimLoadChannels(fileNameAmc, skel);
skelPlayData(skel, channels, frameLength);
