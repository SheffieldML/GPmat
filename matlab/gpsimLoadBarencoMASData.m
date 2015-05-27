function [y, yvar, geneNames, times, scale] = gpsimLoadBarencoMASData

% GPSIMLOADBARENCOMASDATA Load in Martino Barenco's data as processed by MAS5.
% FORMAT

% RETURN y : the normalised expression levels.
% RETURN yvar : the variance of the normalised expression levels.
% RETURN geneNames : the gene names and Affymetrix array tags.
% RETURN times : the times of the expression measurements.
% RETURN scale : the scaling factor applied to normalise.

% SEEALSO : demBarenco2
%
% COPYRIGHT : Pei Gao and Neil D. Lawrence, 2008

% SHEFFIELDML

if exist('./data/barencoMASData.mat') == 2 
  load('./data/barencoMASData.mat');
else
  load('./data/barencoMAS5.mat');
  ind = [2 4 5 1 3];
  for i = 1:5
    geneNames{i} = gene.names{ind(i)};
  end
  
  yFull(1:7,:) = gene.val{1}(:,ind);
  yFull(8:14,:) = gene.val{2}(:,ind);
  yFull(15:21,:) = gene.val{3}(:,ind);
  
  yFullVar(1:7,:) = (gene.var{1}(:,ind)).^2;
  yFullVar(8:14,:) = (gene.var{2}(:,ind)).^2;
  yFullVar(15:21,:) = (gene.var{3}(:,ind)).^2;  
  

  % Rescale so that average standard deviation of curves is 1.
  scale = sqrt(var(yFull));
  scaleMat = (scale'*ones(1,21))';
  yFull = yFull./scaleMat;
  yFullVar = yFullVar./(scaleMat.*scaleMat);
  y{1} = yFull(1:7, :);
  y{2} = yFull(8:14, :);
  y{3} = yFull(15:21, :);
  yvar{1} = yFullVar(1:7, :);
  yvar{2} = yFullVar(8:14, :);
  yvar{3} = yFullVar(15:21, :);
  times = [0 2 4 6 8 10 12]';
  save('./data/barencoMASData.mat', 'y', 'yvar', 'geneNames', 'times', 'scale');
end
