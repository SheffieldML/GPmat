function acclaimPlayFile(fileNameAsf, fileNameAmc, frameLength)

% ACCLAIMPLAYFILE Play motion capture data from a asf and amc file.

% MOCAP

skel = acclaimReadSkel(fileNameAsf);
[channels, skel] = acclaimLoadChannels(fileNameAmc, skel);
skelPlayData(skel, channels, frameLength);
