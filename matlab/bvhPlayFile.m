function bvhPlayFile(fileName, position)

% BVHPLAYFILE Play motion capture data from a bvh file.

% MOCAP

if nargin < 2
  position = 1;
end
[bvhStruct, channels, frameLength] = bvhReadFile(fileName);
bvhPlayData(bvhStruct, channels, frameLength);
