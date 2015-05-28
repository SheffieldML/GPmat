function [y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoData

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

% SHEFFIELDML

if exist('./data/barencoData.mat') == 2 
  load('./data/barencoData.mat');
else
  
  % These excel files include results processed directly from the
  % cel files using the mmgMOS algorithm (Xuejun's code).
  
  % These are the expression levels.
  [numeric1, txt1] = xlsread('./data/resultsMartino_exprs.xls');
  headTxt1 = txt1(1, 2:end);
  tagTxt1 = txt1(2:end, 1);
  
  % These are the standard deviations.
  [numeric2, txt2] = xlsread('./data/resultsMartino_se.xls');
  headTxt2 = txt2(1, 2:end);
  tagTxt2 = txt2(2:end, 1);
  
  if(any(~strcmp(tagTxt2(:), tagTxt1(:))))
    error('Two files are not in same order');
  end
  if(any(~strcmp(headTxt2(:), headTxt1(:))))
    error('Two files are not in same order');
  end
  
  clear gene, clear ind
  % Gene IDs
  % DDB2
  gene{1, 1} = '203409_at';
  gene{1, 2} = 'DDB2';
  % BIK
  gene{2, 1} = '205780_at';
  gene{2, 2} = 'BIK';
  % TNFRSF10b (other tags include 209294_x_at and 210405_x_at)
  gene{3, 1} = '209295_at';
  gene{3, 2} = 'TNFRSF10b';
  % p21 --- we think this is CIp1/p21
  gene{4, 1} = '202284_s_at';
  gene{4, 2} = 'CIp1/p21';
  % p26 --- named as sesn1 in the platform.
  gene{5, 1} = '218346_s_at';
  gene{5, 2} = 'p26 sesn1';
  
  for i = 1:length(gene)
    match = find([strcmp(gene{i, 1}, tagTxt1(:))]);
    if length(match)~=1
      error('Too many or too few matches.');
    else
      ind(i) = match(1);
    end
  end
  order = [1 4 5 6 7 2 3 8 11 12 13 14 9 10 15 18 19 20 21 16 17]; 
  
  % Perform some normalisation.
  % Make sure that the average for each slide in log space is the
  % same.
  mVal = zeros(size(mean(numeric1)));
  mVal = mVal - mean(mVal);
  rawExp = numeric1(ind, order)';
  for i = 1:size(rawExp, 2)
    rawExp(:, i) = rawExp(:, i) - mVal';
  end
  rawVar = numeric2(ind, order)';
  rawVar = rawVar.*rawVar; % convert standard deviations to variances.
  
  yFull = exp(rawExp + rawVar/2);  % Logs are normally distributed
                                   % ... recover mean in exp space.
  yFullVar = (exp(rawVar)-1).*exp(2*rawExp + rawVar); % Logs are
                                                      % normally
                                                      % distributed
                                                      % ... recover
                                                      % variance in exp
                                                      % space.
  

  
  
%     rawExp = zeros(36, length(genes));
%   rawVar = zeros(36, length(genes));
%   yFull = zeros(36, length(genes));
%   yFullVar = zeros(36, length(genes));
%   for k=1:length(genes),
%     I = strcmp(genes{k}, exprs.genes);
%     prof = exprs.data(:, I, :);
%     rawExp(:, k) = squeeze(prof(3, 1, :));
%     rawVar(:, k) = squeeze(diff(prof([4, 2], 1, :)));
%     for l=1:36,
%       t = do_distfit(exp(prof(:, 1, l))', @norminv);
%       yFull(l, k) = t(1);
%       yFullVar(l, k) = t(2) .^ 2;
%     end
%   end
  
  
  
  % Rescale so that average standard deviation of curves is 1.
  scale = mean(sqrt(var(yFull)));
  yFull = yFull/scale;
  yFullVar = yFullVar/(scale*scale);
  y{1} = yFull(1:7, :);
  y{2} = yFull(8:14, :);
  y{3} = yFull(15:21, :);
  yvar{1} = yFullVar(1:7, :);
  yvar{2} = yFullVar(8:14, :);
  yvar{3} = yFullVar(15:21, :);
  times = [0 2 4 6 8 10 12]';
  save('./data/barencoData.mat', 'y', 'yvar', 'gene', 'times', 'scale', 'rawVar', 'rawExp');
end
