function [y, yvar, gene, times, scale, rawExp, rawVar] = gpsimGetDrosData(exprs, exprs_se, genes),

% GPSIMLOADBARENCODATA Load in Martino Barenco's data as processed by mmgMOS.
% FORMAT
% DESC loads in from the two Excel spread sheets
% (resultsMartino_exprs.xls and resultsMartino_se.xls) the data
% from the Barenco et al paper as processed by mmgMOS. 
% RETURN y : the normalised expression levels.
% RETURN yvar : the variance of the normalised expression levels.
% RETURN gene : the gene names and Affymetrix array tags.
% RETURN times : the times of the expression measurements.
% RETURN scale : the scaling factor applied to normalise.
% RETURN rawExp : the raw gene expresion level.
% RETURN rawVar : the raw variance of the gene expression.
% 
% SEEALSO : demBarenco1, demBarencoMap1
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM

if size(genes, 1) < size(genes, 2),
  genes = genes';
end
gene = [genes, genes];

if ~isempty(exprs_se),
  rawExp = zeros(36, length(genes));
  rawVar = zeros(36, length(genes));
  for k=1:length(genes),
    rawExp(:, k) = exprs.(genes{k})';
    rawVar(:, k) = exprs_se.(genes{k})';
  end
  rawVar = rawVar.*rawVar; % convert standard deviations to variances.
  
  yFull = exp(rawExp + rawVar/2);  % Logs are normally distributed
				   % ... recover mean in exp space.
  yFullVar = (exp(rawVar)-1).*exp(2*rawExp + rawVar); % Logs are
						      % normally
						      % distributed
						      % ... recover
						      % variance in exp
						      % space.
else
  rawExp = zeros(36, length(genes));
  rawVar = zeros(36, length(genes));
  yFull = zeros(36, length(genes));
  yFullVar = zeros(36, length(genes));
  for k=1:length(genes),
    I = strcmp(genes{k}, exprs.genes);
    prof = exprs.data(:, I, :);
    rawExp(:, k) = squeeze(prof(3, 1, :));
    rawVar(:, k) = squeeze(diff(prof([4, 2], 1, :)));
    for l=1:36,
      t = do_distfit(exp(prof(:, 1, l))', @norminv);
      yFull(l, k) = t(1);
      yFullVar(l, k) = t(2) .^ 2;
    end
  end
end
  
  
% Rescale so that average standard deviation of curves is 1.
scale = mean(sqrt(var(yFull)));
yFull = yFull/scale;
yFullVar = yFullVar/(scale*scale);
y{1} = yFull(1:12, :);
y{2} = yFull(13:24, :);
y{3} = yFull(25:36, :);
yvar{1} = yFullVar(1:12, :);
yvar{2} = yFullVar(13:24, :);
yvar{3} = yFullVar(25:36, :);
times = (1:12)';
%save('./data/barencoData.mat', 'y', 'yvar', 'gene', 'times', 'scale', 'rawVar', 'rawExp');
