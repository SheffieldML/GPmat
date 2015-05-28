function [y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoPUMAData(option)

% GPSIMLOADBARENCOPUMADATA Load in Martino Barenco's data as re-processed by mmgMOS.
% FORMAT
% DESC loads in from the Excel spread sheets
% the data from the Barenco et al paper as processed by mmgMOS. 
% OPTION defines the way of transforming the data from logged to un-logged space.
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
% COPYRIGHT : Neil D. Lawrence, 2008
% MODIFIED : Pei Gao, 2008
  
% SHEFFIELDML

if nargin < 1
    option = 1;
else
    option = 2;
end

if exist('./data/barencoPUMAData.mat') == 2
    load('./data/barencoPUMAData.mat');
else

    if option == 1
        % These excel files include results processed directly from the
        % cel files using the mmgMOS algorithm (Xuejun's code).

        % These are the expression levels.
%         [numeric1, txt1] = xlsread('./data/barencoPUMA_exprs.xls');
%         headTxt1 = txt1(1, 2:end);
%         tagTxt1 = txt1(2:end, 1);
% 
%         % These are the standard deviations.
%         [numeric2, txt2] = xlsread('./data/barencoPUMA_se.xls');
%         headTxt2 = txt2(1, 2:end);
%         tagTxt2 = txt2(2:end, 1);
        

        % Modified Mauricio Alvarez, to work with the .csv files
        % These are the expression levels.
        s1 = importdata('./data/barencoPUMA_exprs.csv');
        numeric1 =  s1.data;
        txt1 = s1.textdata;
        headTxt1 = txt1(1, 2:end);
        tagTxt1 = txt1(2:end, 1);

        % These are the standard deviations.
        s2 = importdata('./data/barencoPUMA_se.csv');
        numeric2 = s2.data;
        txt2 = s2.textdata;
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

%         for i = 1:length(gene)
%             match = find([strcmp(gene{i, 1}, tagTxt1(:))]);
%             if length(match)~=1
%                 error('Too many or too few matches.');
%             else
%                 ind(i) = match(1);
%             end
%         end
        % Modified Mauricio Alvarez
        for j =1:length(gene),
            ind(j) = strmatch(gene{j,1},tagTxt1, 'exact');
        end
        %order = [1 4 5 6 7 2 3 8 11 12 13 14 9 10 15 18 19 20 21 16 17];
        order = 1:21;
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

        times = [0 2 4 6 8 10 12]';

        yFull = exp(rawExp + rawVar/2);  % Logs are normally distributed
        % ... recover mean in exp space.
        yFullVar = (exp(rawVar)-1).*exp(2*rawExp + rawVar); % Logs are
        % normally
        % distributed
        % ... recover
        % variance in exp
        % space.
    elseif option == 2

        [numeric1, txt1] = xlsread('./data/barencoPUMA_prctile5.xls');
        headTxt1 = txt1(1, 2:end);
        tagTxt1 = txt1(2:end, 1);

        [numeric2, txt2] = xlsread('./data/barencoPUMA_prctile25.xls');
        headTxt2 = txt2(1, 2:end);
        tagTxt2 = txt2(2:end, 1);

        [numeric3, txt3] = xlsread('./data/barencoPUMA_prctile50.xls');
        headTxt3 = txt3(1, 2:end);
        tagTxt3 = txt3(2:end, 1);

        [numeric4, txt4] = xlsread('./data/barencoPUMA_prctile75.xls');
        headTxt4 = txt4(1, 2:end);
        tagTxt4 = txt4(2:end, 1);

        [numeric5, txt5] = xlsread('./data/barencoPUMA_prctile95.xls');
        headTxt5 = txt5(1, 2:end);
        tagTxt5 = txt5(2:end, 1);

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
        %order = [1 4 5 6 7 2 3 8 11 12 13 14 9 10 15 18 19 20 21 16 17];
        order = 1:21;
        % Perform some normalisation.
        % Make sure that the average for each slide in log space is the
        % same.
        mVal = zeros(size(mean(numeric1)));
        mVal = mVal - mean(mVal);
        exprsSet(1,:,:) = numeric1(ind, order);
        exprsSet(2,:,:) = numeric2(ind, order);
        exprsSet(3,:,:) = numeric3(ind, order);
        exprsSet(4,:,:) = numeric4(ind, order);
        exprsSet(5,:,:) = numeric5(ind, order);

        times = [0 2 4 6 8 10 12]';

        for k=1:5
            prof = exprsSet(:, k, :);
            rawExp(:, k) = squeeze(prof(3, 1, :));
            rawVar(:, k) = squeeze(diff(prof([4, 2], 1, :)));
            for l=1:21,
                t = do_distfit(exp(prof(:, 1, l))', @norminv);
                yFull(l, k) = t(1);
                yFullVar(l, k) = t(2) .^ 2;
            end
        end
    end

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

    save('./data/barencoPUMAData.mat', 'y', 'yvar', 'gene', 'times', 'scale', 'rawVar', 'rawExp');
end
