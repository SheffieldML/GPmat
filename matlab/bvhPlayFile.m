function bvhPlayFile(fileName)

% BVHPLAYFILE Play motion capture data from a bvh file.

% MOCAP

[skel, channels, frameLength] = bvhReadFile(fileName);
skelPlayData(skel, channels, frameLength);
