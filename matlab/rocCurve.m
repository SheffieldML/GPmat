function [area, rocPointX, rocPointY] = rocCurve(outputs, labels)

% ROCCURVE Draw ROC curve and return labels.

% NDLUTIL

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