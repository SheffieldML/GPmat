function [area, rocPointX, rocPointY] = rocCurve(outputs, labels)

% ROCCURVE Draw ROC curve and return labels.
%
%	Description:
%
%	[AREA, ROCPOINTX, ROCPOINTY] = ROCCURVE(OUTPUTS, LABELS) draws an
%	ROC curve and returns the area under the ROC curve as well as the
%	points plotted.
%	 Returns:
%	  AREA - teh area under the ROC curve.
%	  ROCPOINTX - the x points of the ROC curve.
%	  ROCPOINTY - the y points of the ROC curve.
%	 Arguments:
%	  OUTPUTS - the outputs from the model (e.g. probabilities of
%	   labels).
%	  LABELS - the true labels associated with the outputs.
%	
%
%	See also
%	TRAPZ


%	Copyright (c) 2004 Neil D. Lawrence



[sOutputs, index] = sort(outputs);
sLabels = labels(index);

for i = length(sOutputs):-1:1;
  % False positives
  rocPointX(length(sOutputs)-i+1) = sum(sLabels(i:end)==-1)/sum(labels==-1);
  rocPointY(length(sOutputs)-i+1) = sum(sLabels(i:end)==1)/sum(labels==1);
end

if nargin < 3
  plot(rocPointX, rocPointY);
end
area = trapz(rocPointX, rocPointY);