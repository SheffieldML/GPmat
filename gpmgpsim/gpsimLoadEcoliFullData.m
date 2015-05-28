function [y, gene, times, scale, rawExp] = gpsimLoadEcoliFullData

% GPSIMLOADECOLIFULLDATA Load in E. coli full data for the repression case.

% RETURN y : the normalised expression levels.
% RETURN yvar : the variance of the normalised expression levels.
% RETURN gene : the gene names and Affymetrix array tags.
% RETURN times : the times of the expression measurements.
% RETURN scale : the scaling factor applied to normalise.
% RETURN rawExp : the raw gene expresion level.
% 
% SEEALSO : demEcoliMap1, gpsimLoadEcoliData
%
% COPYRIGHT : Pei Gao, Neil D. Lawrence and Magnus Rattary, 2008

% SHEFFIELDML

if exist('./data/ecoliData.mat') == 2 
  load('./data/ecoliData.mat');
else
  expData = importdata('./data/ecoliNormalisedData.txt');
  rawExp.data = expData';
  rawExp.genes = {'dinF', 'dinI', 'lexA', 'recA', 'recN', 'ruvA', 'ruvB', ...
                  'sbmC', 'sulA', 'umuC', 'umuD', 'uvrB', 'yebG', 'yjiW' };
  
  targetInd = 1:14;
  
% Perform log-normal transformation.
%   yFull = exp(rawExp.data(:, targetInd));  % Logs are normally distributed
%                                    % ... recover mean in exp space.
  yFull = rawExp.data(:, targetInd); 
  
  % Rescale so that average standard deviation of curves is 1.
  scale = sqrt(var(yFull));
  scaleMat = ones(size(yFull, 1), 1)*scale;
  y{1} = yFull./scaleMat;
  times = [0 5 10 20 40 60]';
  gene = rawExp.genes(targetInd);
%  save('./data/ecoliData.mat', 'y', 'gene', 'times', 'scale', 'rawExp');
end
