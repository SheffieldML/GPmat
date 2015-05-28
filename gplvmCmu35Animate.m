function gplvmCmu35Animate(experimentNo, missingData)

% GPLVMCMU35ANIMATE Animate the test data jointly with predictions.
% FORMAT
% DESC animates the stick figure using 'fill in' from the GP-LVM models.
% ARG experimentNo : the experiment number (1 is four dimensions, 2 is three
% dimensions, 3 is five dimensions). Default value is 1.
% ARG missingData : either 'leg' or 'body' for missing data being the leg
% or the upper body. Default is 'body'.
% 
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : demCmu35gplvm1, demCmu35SequenceOptimisef
  
% FGPLVM
  
dataSetName = 'cmu35gplvm';

% load skeleton
skel = acclaimReadSkel('35.asf');
[tmpchan, skel] = acclaimLoadChannels('35_01.amc', skel);

[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
channelTest = Ytochannels(Ytest);



function channels = Ytochannels(Y)
  
% YTOCHANNELS Convert Y to channel values.

xyzInd = [2];
xyzDiffInd = [1 3];
rotInd = [4 6];
rotDiffInd = [5];
generalInd = [7:38 41:47 49:50 53:59 61:62];
startInd = 1;
endInd = length(generalInd);
channels(:, generalInd) = 180*Y(:, startInd:endInd)/pi;
startInd = endInd + 1;
endInd = endInd + length(xyzDiffInd);
channels(:, xyzDiffInd) = cumsum(Y(:, startInd:endInd), 1);
startInd = endInd + 1;
endInd = endInd + length(xyzInd);
channels(:, xyzInd) = Y(:, startInd:endInd);
startInd = endInd + 1;
endInd = endInd + length(xyzDiffInd);
channels(:, xyzDiffInd) = cumsum(Y(:, startInd:endInd), 1);
startInd = endInd + 1;
endInd = endInd + length(rotInd);
channels(:, rotInd) = asin(Y(:, startInd:endInd))*180/pi;
channels(:, rotInd(end)) = channels(:, rotInd(end))+270;
startInd = endInd + 1;
endInd = endInd + length(rotDiffInd);
channels(:, rotDiffInd) = 0;%cumsum(asin(Y(:, startInd:endInd)), 1))*180/pi;

