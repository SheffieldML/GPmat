function [area, rocPointX, rocPointY] = rocCurve(outputs, labels)

% ROCCURVE Draw ROC curve and return labels.
% FORMAT
% DESC draws an ROC curve and returns the area under the ROC curve
% as well as the points plotted.
% ARG outputs : the outputs from the model (e.g. probabilities of
% labels).
% ARG labels : the true labels associated with the outputs.
% RETURN area : teh area under the ROC curve.
% RETURN rocPointX : the x points of the ROC curve.
% RETURN rocPointY : the y points of the ROC curve.
%
% SEEALSO : trapz
%
% COPYRIGHT : Neil D. Lawrence, 2004


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
