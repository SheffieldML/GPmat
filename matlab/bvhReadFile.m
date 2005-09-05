function [bvhStruct, channels, frameLength] = bvhReadFile(fileName)

% BVHREADFILE Reads a bvh file into a tree structure.

% MOCAP

% a regular expression for floats
numPat = '(-?[0-9]*\.?[0-9]*)';
% a regular expression for positive ints
intPat = '([0-9]+)';


fid = fopen(fileName, 'r');
numBracket = 0;
% start the backtrace with 0 (the parent of the root node).
backTrace = 0;
jointNo = 0;
channelNo = 0;
while ~feof(fid)
  lin = fgetl(fid);
  [token, R] = strtok(lin);

  switch token
   
   case {'ROOT', 'JOINT', 'End'}
    % A joint is starting
    parentNo = backTrace(numBracket+1);
    jointNo = jointNo + 1;
    name = strtok(R); 
    bvhStruct(jointNo) = struct('name', name, ...
                                'offset', [], ...
                                'channels', [], ...
                                'parent', parentNo, ...
                                'order', [], ...
                                'rotInd', [], ...
                                'posInd', [], ...
                                'children', []);
    backTrace = [backTrace jointNo];
    
   case 'OFFSET'
    % Read in the offset.    
    offsets = regexp(R, numPat, 'tokens');
    if length(offsets) == 3
      for i = 1:3
        offset = str2num(offsets{i}{1});
        bvhStruct(jointNo).offset(1, i) = offset;
      end
    else 
      error('Offset is incorrect length.');
    end
   
   case '{'    
    numBracket = numBracket + 1;
    
   case '}'
    numBracket = numBracket - 1;
    backTrace = backTrace(1:end-1);
   
   case 'CHANNELS'
    % Read in the channels     
    [tok, R] = strtok(R); 
    numChannels = str2num(tok);
    orderNo = 0;
    for i = 1:numChannels 
      [tok, R] = strtok(R);
      bvhStruct(jointNo).channels{i} = tok;
      switch tok
       case 'Xrotation'
        orderNo = orderNo + 1;
        bvhStruct(jointNo).rotInd(1, 1) = channelNo + i;
        bvhStruct(jointNo).order(1, orderNo) = 'x';
       case 'Yrotation'
        orderNo = orderNo + 1;
        bvhStruct(jointNo).rotInd(1, 2) = channelNo + i;
        bvhStruct(jointNo).order(1, orderNo) = 'y';
       case 'Zrotation'
        orderNo = orderNo + 1;
        bvhStruct(jointNo).rotInd(1, 3) = channelNo + i;
        bvhStruct(jointNo).order(1, orderNo) = 'z';
       case 'Xposition'
        bvhStruct(jointNo).posInd(1, 1) = channelNo + i;
       case 'Yposition'
        bvhStruct(jointNo).posInd(1, 2) = channelNo + i;
       case 'Zposition'
        bvhStruct(jointNo).posInd(1, 3) = channelNo + i;
      end        
    end
    bvhStruct(jointNo).order = char(bvhStruct(jointNo).order);
    channelNo = channelNo + numChannels;
   
   case 'MOTION'
    lin = fgetl(fid);
    pat = ['Frames:\s*' intPat]; 
    frames = regexpi(lin, pat, 'tokens');
    if length(frames) == 1
      dataStruct.numFrames = str2num(frames{1}{1});
    else
      error('Cannot determine number of frames in file.')
    end
    
    lin = fgetl(fid);
    pat = ['Frame Time:\s*' numPat]; 
    frameLength = regexpi(lin, pat, 'tokens');
    if length(frameLength)
      frameLength = str2num(frameLength{1}{1});
    else
      error('Cannot determine length of frames in file.')
    end
    
    channels = [];
    while ~feof(fid)
      lin = fgetl(fid);
      chanLine = sscanf(lin, '%f');
      channels = [channels; chanLine'];
    end
  end
end
fclose(fid);
bvhStruct = findChildren(bvhStruct);

channels = channelsAngles(channels, bvhStruct);

function tree = findChildren(tree)

% FINDCHILDREN Add the children to the tree structure.

for i = 1:length(tree)
  for j = 1:length(tree(i).parent)
    if tree(i).parent(j)
      tree(tree(i).parent(j)).children ...
          = [tree(tree(i).parent(j)).children i];
    end
  end
end


function channels = channelsAngles(channels, bvhStruct);

% CHANNELSANGLES Try and remove artificial discontinuities associated with angles.

for i=1:length(bvhStruct)
  for j=1:length(bvhStruct(i).rotInd)
    col = bvhStruct(i).rotInd(j);
    for k=2:size(channels, 1)
      diff=channels(k, col)-channels(k-1, col);
      if abs(diff+360)<abs(diff)
        channels(k:end, col)=channels(k:end, col)+360;
      elseif abs(diff-360)<abs(diff)
        channels(k:end, col)=channels(k:end, col)-360;
      end
    end
  end
end