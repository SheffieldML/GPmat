function [channels, skel] = acclaimSemiLoadChannels(fileName, skel, iSubset)

% ACCLAIMSEMILOADCHANNELS Load the channels from an AMC file for a subset
% of the frames.
% FORMAT
% DESC loads channels from an acclaim motion capture file.
% ARG fileName : the file name to load in.
% ARG skel : the skeleton structure for the file.
% RETURN channels : the channels read in for the file.
% RETURN skel : the skeleton for the file.
% 
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Alfredo A. Kalaitzis, 2012
% 
% SEEALSO : acclaimLoadChannels, acclaimReadSkel

% MOCAP

if nargin < 3
    iSubset = 1:acclaimNumberOfFrames(fileName);
elseif iscell(iSubset)
    iSubset = cell2mat(iSubset);
end

fid=fopen(fileName, 'rt');

lin = getline(fid);
lin = strtrim(lin);
while ~strcmp(lin,':DEGREES')
  lin = getline(fid);
  lin = strtrim(lin);
end

counter = 0;

if ~isempty(iSubset)
    while ~strcmp(lin, num2str(iSubset(1)))
        lin = getline(fid);
    end
    iSubset = iSubset(2:end);
end

while lin~=-1
  lin = strtrim(lin);
  parts = tokenise(lin, ' ');
  if length(parts)==1
    frameNo = str2num(parts{1});
    if ~isempty(frameNo)
      counter = counter + 1;
      if counter ~= frameNo
        error('Unexpected frame number.');
      end
    else
      error('Single bone name  ...');
    end
  else
    ind = skelReverseLookup(skel, parts{1});
    for i = 1:length(skel.tree(ind).channels)
      bone{ind}{frameNo}(i) = str2num(parts{i+1});
    end
  end
  lin = getline(fid);
end
fclose(fid);
numFrames = counter;
numChannels = 0;
for i = 1:length(skel.tree)
  numChannels = numChannels + length(skel.tree(i).channels);
end
channels = zeros(numFrames, numChannels);

endVal = 0;
for i =1:length(skel.tree)
  if length(skel.tree(i).channels)>0
    startVal = endVal + 1;
    endVal = endVal + length(skel.tree(i).channels);
    channels(:, startVal:endVal)= reshape([bone{i}{:}], ...
                                          length(skel.tree(i).channels), numFrames)';
    
  end
  [skel.tree(i).rotInd, skel.tree(i).posInd] = ...
      resolveIndices(skel.tree(i).channels, startVal);
end

channels = smoothAngleChannels(channels, skel);

function [rotInd, posInd] = resolveIndices(channels, startVal)

% RESOLVEINDICES Get indices from the channels.

baseChannel = startVal - 1;
rotInd = zeros(1, 3);
posInd = zeros(1, 3);
for i = 1:length(channels)
  switch channels{i}
   case 'Xrotation'
    rotInd(1, 1) = baseChannel + i;
   case 'Yrotation'
    rotInd(1, 2) = baseChannel + i;
   case 'Zrotation'
    rotInd(1, 3) = baseChannel + i;
   case 'Xposition'
    posInd(1, 1) = baseChannel + i;
   case 'Yposition'
    posInd(1, 2) = baseChannel + i;
   case 'Zposition'
    posInd(1, 3) = baseChannel + i;
  end        
end



