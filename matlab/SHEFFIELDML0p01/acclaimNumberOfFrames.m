function nFrames = acclaimNumberOfFrames(fileName)

% ACCLAIMNUMBEROFFRAMES Extract the number of frames.
%
%	Description:
%
%	NFRAMES = ACCLAIMNUMBEROFFRAMES(FILENAME) counts the total number of
%	frames from an acclaim motion capture file.
%	 Returns:
%	  NFRAMES - the total number of frames in the file.
%	 Arguments:
%	  FILENAME - the file name to load in.
%	
%
%	See also
%	


%	Copyright (c) 2012 Alfredo A. Kalaitzis


fid = fopen(fileName, 'rt');

% Assumes the number of lines in each frame is fixed.
a = cell(textscan(fid, '%s'));
a = a{1};
nFrames = str2double(a{end-91}); % 91 strings from eof to last frame #

fclose(fid);