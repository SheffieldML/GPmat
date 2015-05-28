function [X, y, XTest, yTest] = mapLoadData(dataset, seedVal)

% MAPLOADDATA Load a mapping model dataset (e.g. classification, regression).
% FORMAT
% DESC loads a data set for a mapping modelling problem
% (e.g. classification or regression).
% ARG dataset : the name of the data set to be loaded. Currently
% the possible names are 'usps', 'usps.' where . is a digit between
% 0 and 9, 'spgp1d', 'twoclusters', 'regressionOne',
% 'regressionTwo', 'regressionThree', 'regressionFour',
% 'classificationTwo', 'classificationThree', 'orderedOne',
% 'orderedTwo', 'ionosphere'.
% ARG seedVal : a seed value for generating the data set (default
% is 1e5). Note that in many cases a generated data set will be
% saved and loaded from disc the next time the function is
% called. In these cases this value, if changed, won't have an effect.
% RETURN X : the training input data loaded in.
% RETURN y : the training target data loaded in.
% RETURN XTest : the test input data loaded in. If no test set is
% available it is empty.
% RETURN yTest : a test target data.
%
% SEEALSO : lvmLoadData, datasetsDirectory
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% DATASETS

if nargin < 2
    seedVal = 1e5;
end
randn('seed', seedVal)
rand('seed', seedVal)
XTest = [];
yTest = [];

baseDir = datasetsDirectory;
dirSep = filesep;
semiSup = false;

if length(dataset)>5 && strcmp(dataset(1:6), 'gunnar')
    % Data set is one of Gunnar Raetsch's
    ind = find(dataset==':');
    dataSetName = dataset(ind(1)+1:ind(2)-1);
    dataSetNum = dataset(ind(2)+1:end);
    filebase = [baseDir dirSep 'gunnar' dirSep dataSetName dirSep dataSetName];
    X=load([filebase '_train_data_' num2str(dataSetNum) '.asc']);
    y=load([filebase '_train_labels_' num2str(dataSetNum) '.asc']);
    XTest=load([filebase '_test_data_' num2str(dataSetNum) '.asc']);
    yTest=load([filebase '_test_labels_' num2str(dataSetNum) '.asc']);
    return
elseif length(dataset)>3 && strcmp(dataset(1:4), 'semi')
    % Data set is semi-supervised learning.
    ind = find(dataset==':');
    labProb = str2num(dataset(ind(2)+1:end));
    % Extract dataset part
    dataset = dataset(ind(1)+1:ind(2)-1);
    semiSup = true;
end

switch dataset

    case 'silhouette'
        % Ankur Agarwal and Bill Trigg's silhoutte data.
        load([baseDir dirSep 'mocap' dirSep 'ankur' dirSep 'ankurDataPoseSilhouette']);
        inMean = mean(Y);
        inScales = sqrt(var(Y));
        X = Y - repmat(inMean, size(Y, 1), 1);
        X = X./repmat(inScales, size(Y, 1), 1);

        XTest = Y_test - repmat(inMean, size(Y_test, 1), 1);
        XTest = XTest./repmat(inScales, size(Y_test, 1), 1);
        y = Z;
        yTest = Z_test;


    case 'cedar69'
        % Data for ICML 2001 paper on noisy KFD.
        digOne = 6;
        digTwo = 9;

        directory = [baseDir '\cedar.cd\matlab\16x16\'];
        load([directory 'digit' num2str(digOne)]);
        X = digitData;
        numDataOne = size(digitData, 1);
        y = zeros(size(digitData, 1), 1);
        mydigitData = [];
        digitData = [];
        load([directory 'digit' num2str(digTwo)]);
        X = [X; digitData];
        X = double(X);
        X = -(X-128)/128;
        numDataTwo = size(digitData, 1);
        y = [y; ones(size(digitData, 1), 1)];

        if nargin>2
            % Load test data
            load([directory 'bsdigit' num2str(digOne)]);
            XTest = digitData;
            numDataOne = size(digitData, 1);
            yTest = zeros(size(digitData, 1), 1);
            digitData = [];
            load([directory 'bsdigit' num2str(digTwo)]);
            XTest = [XTest; digitData];
            XTest = double(XTest);
            XTest = -(XTest-128)/128;
            yTest = [yTest; ones(size(digitData, 1), 1)];
        end

    case 'sky'
        load([baseDir 'skydatatrain.mat']);
        X = double(x)/255;
        y = double(t)*2 -1;
        load([baseDir 'skydatavalid.mat']);
        XTest = double(x)/255;
        yTest = double(t)*2 -1;

    case 'unlabelledOne'
        numDataPart = 100;
        [X, y] = generateCrescentData(numDataPart);

    case 'classificationOne'
        try
            load([baseDir 'classificationOneData'])
        catch
            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
                y = [ones(200, 1); -ones(100, 1)];
                save([baseDir 'classificationOneData.mat'], 'X', 'y')
            else
                error(lasterr);
            end
        end
    case 'usps'
        load([baseDir 'usps_train']);
        X = ALL_DATA;
        range =  min(ALL_T):max(ALL_T);
        for i = 1:length(range)
            y(:, i) = (ALL_T == range(i))*2 - 1;
        end
        if nargout > 2
            load([baseDir 'usps_test']);
            XTest = ALL_DATA;
            range =  min(ALL_T):max(ALL_T);
            for i = 1:length(range)
                yTest(:, i) = (ALL_T == range(i))*2 - 1;
            end
        end

    case {'usps0', 'usps1', 'usps2', 'usps3', 'usps4', 'usps5', 'usps6', 'usps7', 'usps8', 'usps9'}
        digitNo = str2num(dataset(end));
        load([baseDir 'usps_train'])
        X = ALL_DATA;
        range =  min(ALL_T):max(ALL_T);
        for i = 1:length(range)
            y(:, i) = (ALL_T == range(i))*2 - 1;
        end
        if nargout > 2
            load([baseDir 'usps_test']);
            XTest = ALL_DATA;
            range =  min(ALL_T):max(ALL_T);
            for i = 1:length(range)
                yTest(:, i) = (ALL_T == range(i))*2 - 1;
            end
        end
        y = y(:, digitNo+1);
        yTest = yTest(:, digitNo+1);
    case 'threeFive'
        load([baseDir 'usps_train'])
        X = ALL_DATA;
        y = ALL_T;
        load([baseDir 'usps_test'])
        XTest = ALL_DATA;
        yTest = ALL_T;
        classTrue = 3;
        for i = [0 1 2 4 6 7 8 9];
            index = find(y == i);
            X(index, :) = [];
            y(index, :) = [];
            index = find(yTest == i);
            XTest(index, :) = [];
            yTest(index, :) = [];
        end
        y = (y == classTrue)*2 - 1;
        yTest = (yTest == classTrue)*2 - 1;
    case 'fourNine'
        load([baseDir 'usps_train'])
        X = ALL_DATA;
        y = ALL_T;
        load([baseDir 'usps_test'])
        XTest = ALL_DATA;
        yTest = ALL_T;
        classTrue = 4;
        for i = [0 1 2 3 5 6 7 8];
            index = find(y == i);
            X(index, :) = [];
            y(index, :) = [];
            index = find(yTest == i);
            XTest(index, :) = [];
            yTest(index, :) = [];
        end
        y = (y == classTrue)*2 - 1;
        yTest = (yTest == classTrue)*2 - 1;

    case 'thorsten'
        [y, X] = svmlread([baseDir 'example2/train_transduction.dat']);
        [yTest, XTest] = svmlread([baseDir 'example2/test.dat']);

    case 'spgp1d'
        try
            load([baseDir 'spgp1DData.mat'])
        catch
            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 1;
                N = 500;
                X = 2*rand(N, numIn)-1;
                kern = kernCreate(X, 'rbf');
                kern.variance = 1;
                kern.inverseWidth = 20;
                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.1;

                save([baseDir 'spgp1DData.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr);
            end
        end
    case 'twoclusters'
        try
            load([baseDir 'twoclusters.mat'])
        catch
            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 1;
                N = 200;
                X1 = rand(N/2, numIn)-1.5;
                X2 = rand(N/2, numIn)+.5;
                X = [X1; X2];
                kern = kernCreate(X, 'rbf');
                kern.variance = 1;
                kern.inverseWidth = 100;
                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.1;

                save([baseDir 'twoclusters.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr);
            end
        end
    case 'regressionOne'
        try
            load([baseDir 'regressionOneData.mat'])
        catch

            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;
                X = zeros(N, numIn);
                X(1:floor(N/2), :) = ...
                    randn(floor(N/2), numIn)*.5+1;
                X(floor(N/2)+1:end, :) = ...
                    randn(ceil(N/2), numIn)*.5-1;
                kern = kernCreate(X, 'rbfard');
                kern.variance = 1;
                kern.inverseWidth = 20;
                kern.inputScales = [0 0.999];

                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;

                save([baseDir 'regressionOneData.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr);
            end

        end

    case 'regressionTwo'
        try
            load([baseDir 'regressionTwoData.mat'])
        catch

            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;
                X = zeros(N, numIn);
                X(1:floor(N/2), :) = ...
                    randn(floor(N/2), numIn)*.5+1;
                X(floor(N/2)+1:end, :) = ...
                    randn(ceil(N/2), numIn)*.5-1;
                kern = kernCreate(X, 'rbfard');
                kern.variance = 1;
                kern.inverseWidth = 20;
                kern.inputScales = [0.999 .2];

                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;

                save([baseDir 'regressionTwoData.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr)
            end

        end

    case 'regressionThree'
        try
            load([baseDir 'regressionThreeData.mat'])
        catch

            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;
                X = zeros(N, numIn);
                X(1:floor(N/2), :) = ...
                    randn(floor(N/2), numIn)*.5+1;
                X(floor(N/2)+1:end, :) = ...
                    randn(ceil(N/2), numIn)*.5-1;
                kern = kernCreate(X, 'lin');
                kern.variance = 1;

                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
                save([baseDir 'regressionThreeData.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr);
            end
        end

    case 'regressionFour'
        try
            load([baseDir 'regressionFourData.mat'])
        catch

            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;
                X = zeros(N, numIn);
                X(1:floor(N/2), :) = ...
                    randn(floor(N/2), numIn)*.5+1;
                X(floor(N/2)+1:end, :) = ...
                    randn(ceil(N/2), numIn)*.5-1;
                kern = kernCreate(X, 'mlp');
                kern.variance = 1;
                kern.weightVariance = 1;
                kern.biasVariance = 1;
                K = kernCompute(kern, X);
                y = real(gaussSamp(K, 1)') + randn(N, 1)*0.01;
                save([baseDir 'regressionFourData.mat'], 'numIn', 'N', 'X', 'y')
            else
                error(lasterr);
            end
        end

    case 'classificationTwo'

        try
            load([baseDir 'classificationTwoData.mat'])
        catch
            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;

                X = zeros(N, numIn);
                X = rand(N, numIn);

                kern = kernCreate(X, 'rbf');
                kern.variance = 10;
                kern.inverseWidth = 10;

                K = kernCompute(kern, X);
                u = real(gaussSamp(K, 1)');

                p = cumGaussian(u);
                y = 2*(rand(size(u))>p)-1;
                save([baseDir 'classificationTwoData.mat'], 'numIn', 'N', 'X', 'u', 'y')
            else
                error(lasterr);
            end

        end

    case 'classificationThree'

        try
            load([baseDir 'classificationThreeData.mat'])
        catch
            [void, errid] = lasterr;
            if strcmp(errid, 'MATLAB:load:couldNotReadFile')
                numIn = 2;
                N = 500;

                X = zeros(N, numIn);
                X = rand(N, numIn);

                kern = kernCreate(X, 'rbf');
                kern.variance = 10;
                kern.inverseWidth = 10;

                K = kernCompute(kern, X);
                u = real(gaussSamp(K, 1)');
                a = 3;
                pMinus = cumGaussian(u-a/2);
                pPlus = cumGaussian(-u-a/2);
                p =rand(size(u));
                indMinus = find(p<pMinus);
                indPlus = find(p>pMinus && p<pMinus+pPlus);
                indNone = find(p>pMinus+pPlus);
                y = zeros(N, 1);
                y(indPlus) = 1;
                y(indMinus) = -1;
                y(indNone, :) = [];
                X(indNone, :) = [];
                save([baseDir 'classificationThreeData.mat'], 'numIn', 'N', 'X', 'u', 'y')
            else
                error(lasterr);
            end

        end

    case 'orderedOne'
        dataPerCat = 30;
        spacing = 3;

        % Generate a toy data-set of (linear) ordered categories.
        X = [randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat, 2); ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]];
        y = [zeros(dataPerCat, 1); ...
            repmat(1, dataPerCat, 1); repmat(2, dataPerCat, 1); ...
            repmat(3, dataPerCat, 1); repmat(4, dataPerCat, 1); ...
            repmat(5, dataPerCat, 1); repmat(6, dataPerCat, 1)];

    case 'orderedTwo'
        dataPerCat = 30;
        spacing = 3;

        % Generate a toy data-set of (linear) ordered categories.
        thetaR = [randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)-[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat, 2); ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(2*spacing, dataPerCat, 1)]; ...
            randn(dataPerCat,2)+[zeros(dataPerCat, 1) repmat(3*spacing, dataPerCat, 1)]];
        thetaR(:, 1) = thetaR(:, 1);
        thetaR(:, 2) = thetaR(:, 2) + 15;
        X = [sin(thetaR(:, 1)).*thetaR(:, 2) cos(thetaR(:, 1)).*thetaR(:, 2)];
        y = [zeros(dataPerCat, 1); ...
            repmat(1, dataPerCat, 1); repmat(2, dataPerCat, 1); ...
            repmat(3, dataPerCat, 1); repmat(4, dataPerCat, 1); ...
            repmat(5, dataPerCat, 1); repmat(6, dataPerCat, 1)];


    case 'ionosphere'
        X = zeros(351, 34);
        fid = fopen([baseDir 'ionosphere.data'], 'r');
        lin = getline(fid);
        i = 0;
        while(lin ~= -1)
            i = i+1;
            elements = tokenise(lin, ',');
            for j = 1:length(elements)-1
                X(i, j) = str2num(elements{j});
            end
            switch(elements{end})
                case 'g'
                    y(i, 1) = 1;
                case 'b'
                    y(i, 1) = -1;
            end
            lin = getline(fid);
        end
        ind = randperm(351);
        trainInd = ind(1:200);
        testInd = ind(201:end);
        XTest = X(testInd, :);
        yTest = y(testInd, :);
        X = X(trainInd, :);
        y = y(trainInd, :);

        %/~
    case 'pumadynSeeger'
        data = load([baseDir 'Dataset.data']);
        ind = randperm(size(data, 1));
        indTr = ind(1:7168);
        indTe = ind(7169:end);
        X = data(indTr, 1:end-1);
        y = data(indTr, end);
        XTest = data(indTe, 1:end-1);
        yTest = data(indTe, end);
        Xscale = sqrt(var(X));
        for j = 1:length(Xscale);
            X(:, j) = X(:, j)/Xscale(j);
            XTest(:, j) = XTest(:, j)/Xscale(j);
        end
        augX = [X ones(size(X, 1), 1)];
        w = inv(augX'*augX)*augX'*y;
        y = y -augX*w;
        augTestX = [XTest ones(size(XTest, 1), 1)];
        yTest = yTest - augTestX*w;
        yscale = sqrt(var(y));
        y = y/yscale;
        yTest = yTest/yscale;
    case 'pumadyn'

        % Data is variance 1, no need to normalise.
        data = load([baseDir 'Dataset.data']);
        ind = randperm(size(data, 1));
        indTr = ind(1:7168);
        indTe = ind(7169:end);
        X = data(indTr, 1:end-1);
        y = data(indTr, end);
        XTest = data(indTe, 1:end-1);
        yTest = data(indTe, end);
        %~/
        %/~
    case 'lfmOde'
        data = load([baseDir 'datasetODE31_150_5Ind0']);
        X = data.Xtrain{1};
        y = cell2mat(data.Ytrain');
        XTest = data.Xtest{1};
        yTest = cell2mat(data.Ytest');

    case 'lfmOdeInd'
        data = load([baseDir 'datasetODE31_100_5Ind1']);
        X = data.Xtrain{1};
        y = cell2mat(data.Ytrain');
        XTest = data.Xtest{1};
        yTest = cell2mat(data.Ytest');
        %~/
    case 'ggToy'
        try
            load([baseDir 'toyMultigp1D.mat']);
        catch
            randn('seed', 1e5); % This was taken to make the outputs equal to the ones of the NIPS paper
            rand('seed', 1e5);  % This was taken to make the outputs equal to the ones of the NIPS paper
            nout = 4;
            nlf = 1;
            d = nout + nlf;
            precisionG = [50 50 300 200];
            variance = [1 1 5 5];
            precisionU = 100;
            sigma2Latent = 1;
            N = 500;
            N2 = 100;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            XX = cell(1, d);
            XX{1} = x2;
            for i =1: nout;
                XX{nlf+i} = x;
            end
            kernType{1} = multigpKernComposer('gg', d, nlf, 'ftc', 1);
            kern = kernCreate(XX,  kernType{:});
            kern.comp{1}.precisionU = precisionU;
            kern.comp{1}.sigma2Latent = sigma2Latent;
            for k = 1:nout,
                kern.comp{1+k}.precisionU = precisionU;
                kern.comp{1+k}.sigma2Latent = sigma2Latent;
                kern.comp{1+k}.precisionG = precisionG(k);
                kern.comp{1+k}.sensitivity = variance(k);
            end
            K = kernCompute(kern, XX);
            yu = gsamp(zeros(size(K,1),1), K, 1);
            y = yu(size(XX{1},1)+1:end);
            Y = reshape(y,size(XX{2},1),nout);
            F = reshape(y,size(XX{2},1),nout);
            for k=1:nout,
                F(:,k) = Y(:,k);
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            X = cell(1,nout);
            y = cell(1,nout);
            f = cell(1,nout);
            for k =1:nout,
                X{k} = XX{k+nlf};
                y{k} = Y(:,k);
                f{k} = F(:,k);
            end
            XTest = [];
            yTest = f;
            save([baseDir 'toyMultigp1D.mat'], 'X', 'y', 'XTest', 'yTest')
        end

    case 'ggToyTrainTest'
        try
            load([baseDir 'toyMultigp1DTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        catch
            [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData('ggToy');
            nout = size(XTemp,2);
            ntrainx = 200;
            randn('seed', 1e4); % This was taken to make the outputs equal to the ones of the NIPS paper
            rand('seed', 1e4);  % This was taken to make the outputs equal to the ones of the NIPS paper
            maxl = length(XTemp{1});
            X = cell(1,nout);
            y = cell(1,nout);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = XTemp{k}(pindx,:);
                y{k} = yTemp{k}(pindx,:);
                yTest{k} = yTemp{k}(indx(ntrainx+1:end),:);
            end
            XTest = XTemp{1}(indx(ntrainx+1:end),:)';
            yTest{end+1} = yTestTemp; %It has the groundtruth
            save([baseDir 'toyMultigp1DTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        end



    case 'ggToyMissing'
        try
            load([baseDir 'toyMultigp1DMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        catch
            [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData('ggToy');
            nout = size(XTemp,2);
            missingData = cell(nout,1);
            missingData{1} = 101:160;
            missingData{4} = 21:90;
            ntrainx = 200;
            randn('seed', 1e4); % This was taken to make the outputs equal to the ones of the NIPS paper
            rand('seed', 1e4);  % This was taken to make the outputs equal to the ones of the NIPS paper
            maxl = length(XTemp{1});
            X = cell(1,nout);
            y = cell(1,nout);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = XTemp{k}(pindx,:);
                X{k}(missingData{k},:)= [];
                y{k} = yTemp{k}(pindx,:);
                y{k}(missingData{k}) = [];
                yTest{k} = yTemp{k}(indx(ntrainx+1:end),:);
            end
            XTest = XTemp{1}(indx(ntrainx+1:end),:)';
            yTest{end+1} = yTestTemp; % It has the groundtruth
            save([baseDir 'toyMultigp1DMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToy'
        try
            load([baseDir 'toywhiteMultigp1D.mat']);
        catch
            randn('seed', 1e5); % This was taken to make the outputs equal to the ones of the NIPS paper
            rand('seed', 1e5);  % This was taken to make the outputs equal to the ones of the NIPS paper
            nout =4;
            nlf = 1;
            d = nlf + nout;
            precisionG = [50 50 300 200];
            variance = [1 1 5 5];
            N = 500;
            x = linspace(-1,1,N)';
            kernType{1} = multigpKernComposer('ggwhite', d, nlf, 'ftc', 1);
            kern = kernCreate(x,  kernType{:});
            kern.comp{1}.variance = 1;
            for k = 1:nout,
                kern.comp{1+k}.precisionG = precisionG(k);
                kern.comp{1+k}.variance = variance(k);
            end
            K = kernCompute(kern, x);
            yu = gsamp(zeros(size(K,1),1), K, 1);
            y = yu(size(x,1)+1:end);
            Y = reshape(y,size(x,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            X = cell(1,nout);
            y = cell(1,nout);
            for k =1:nout,
                X{k} = x;
                y{k} = Y(:,k);
            end
            XTest = [];
            yTest = [];
            save([baseDir 'toywhiteMultigp1D.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToyMissing'
        try
            load([baseDir 'toywhiteMultigp1DMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        catch
            [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData('ggwhiteToy');
            nout = size(XTemp,2);
            missingData = cell(nout,1);
            missingData{1} = 101:160;
            missingData{4} = 21:90;
            ntrainx = 200;
            randn('seed', 1e4); % This was taken to make the outputs equal to the ones of the NIPS paper
            rand('seed', 1e4);  % This was taken to make the outputs equal to the ones of the NIPS paper
            maxl = length(XTemp{1});
            X = cell(1,nout);
            y = cell(1,nout);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = XTemp{k}(pindx,:);
                X{k}(missingData{k},:)= [];
                y{k} = yTemp{k}(pindx,:);
                y{k}(missingData{k}) = [];
                yTest{k} = yTemp{k}(indx(ntrainx+1:end),:);
            end
            XTest = XTemp{1}(indx(ntrainx+1:end),:)';
            save([baseDir 'toywhiteMultigp1DMissing.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToyHighDim'
        try
            load([baseDir 'ggwhiteToyHighDim.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 1;
            inputDim = 2;
            prec = [50 10 5 100 30];
            sens_y_noise = [1 2 2.2 1.5 1];
            N = 20;
            N2 = 20;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            [Xout1, Xout2 ] = ndgrid(x);
            [Xlf1, Xlf2 ] = ndgrid(x2);
            Xout = [Xout1(:) Xout2(:)];
            Xlf = [Xlf1(:)  Xlf2(:)];
            Kyy = cell(nout);
            ggKern1Noise.inputDimension = inputDim;
            ggKern1Noise.sigma2Noise = noise;
            ggKern2Noise.inputDimension = inputDim;
            for i = 1:nout,
                ggKern1Noise.precisionG = prec(i)*ones(inputDim,1);
                ggKern1Noise.variance = sens_y_noise(i);
                for j = 1:nout,
                    ggKern2Noise.precisionG = prec(j)*ones(inputDim,1);
                    ggKern2Noise.variance = sens_y_noise(j);
                    Kyy{i,j} = ggwhiteXggwhiteKernCompute(ggKern1Noise, ggKern2Noise, Xout);
                end
            end
            Kyu = cell(nout, nin);
            ggKernNoise.inputDimension = inputDim;
            ggKernNoise.variance = noise;
            for i = 1:nout,
                ggKern1Noise.variance = sens_y_noise(i);
                ggKern1Noise.precisionG = prec(i)*ones(inputDim,1);
                Kyu{i,1} = ggwhiteXwhiteKernCompute(ggKern1Noise, ggKernNoise, Xout, Xlf);
            end
            Kuu = cell(nin,1);
            Kuu{1} = noise*eye(size(Xlf,1));
            K = [cell2mat(Kuu) cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = gsamp(zeros(size(K,1),1), K, 1);
            u = yu(1:size(Xlf,1));
            y = yu(size(Xlf,1)+1:end);
            U = reshape(u,size(Xlf,1),nin);
            Y = reshape(y,size(Xout,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =200;
            maxl = size(Xout,1);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = Xout(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = Xout(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            X{nout+1} = Xlf;
            save([baseDir 'ggwhiteToyHighDim.mat'], 'X', 'y', 'XTest', 'yTest')
        end

    case 'ggwhiteToyMissingHighDim'
        try
            load([baseDir 'ggwhiteToyMissingHighDim.mat']);
        catch
            [X, y, XTest, yTest] = mapLoadData('ggwhiteToyHighDim');
            nout = 5;
            missingData = cell(nout,1);
            missingData{1} = 46:108;
            missingData{4} = 147:184;
            for k =1:nout,
                X{k}(missingData{k},:)= [];
                y{k}(missingData{k}) = [];
            end
            save([baseDir 'ggwhiteToyMissingHighDim.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'simToyCombined'
        try
            load([baseDir 'simToyCombined.mat']);
        catch
            noise = 1;
            nout =5;
            nin = 2;
            decay = [2 10 5 0.1 7.5];
            variance = [1.5 1.5 4 10 2];
            sensitivity = [1 0.3 1.2 0.2 0.5];
            inverseWidth = 100;
            N = 200;
            N2 = 50;
            x = linspace(0,1,N)';
            x2 = linspace(0,1,N2)';
            Kyy = cell(nout);
            simKern1.inputDimension = size(x,2);
            simKern1 = simKernParamInit(simKern1);
            simKern1.inverseWidth = inverseWidth;
            simwhiteKern1.inputDimension =size(x,2);
            simwhiteKern1 = simwhiteKernParamInit(simwhiteKern1);
            simwhiteKern1.variance = noise;
            simKern2.inputDimension = size(x,2);
            simKern2 = simKernParamInit(simKern2);
            simKern2.inverseWidth = inverseWidth;
            simwhiteKern2.inputDimension = size(x,2);
            simwhiteKern2 = simwhiteKernParamInit(simwhiteKern2);
            simwhiteKern2.variance = noise;
            for i = 1:nout,
                simKern1.variance = variance(i);
                simKern1.decay = decay(i);
                simwhiteKern1.decay = decay(i);
                simwhiteKern1.sensitivity = sensitivity(i);
                for j = 1:nout,
                    simKern2.variance = variance(j);
                    simKern2.decay = decay(j);
                    simwhiteKern2.decay = decay(j);
                    simwhiteKern2.sensitivity = sensitivity(j);
                    Kyy{i,j} = simXsimKernCompute(simKern1, simKern2, x) + ...
                        simwhiteXsimwhiteKernCompute(simwhiteKern1, simwhiteKern2, x);
                end
            end
            Kyu = cell(nout, 2);
            rbfKern.inputDimension = size(x,2);
            rbfKern = rbfKernParamInit(rbfKern);
            rbfKern.inverseWidth = inverseWidth;
            for i = 1:nout,
                simKern1.variance = variance(i);
                simKern1.decay = decay(i);
                Kyu{i,1} = simXrbfKernCompute(simKern1, rbfKern, x, x2);
            end
            whiteKern.inputDimension = size(x,2);
            whiteKern = whiteKernParamInit(whiteKern);
            whiteKern.variance = noise;
            for i = 1:nout,
                simwhiteKern1.decay = decay(i);
                simwhiteKern1.sensitivity = sensitivity(i);
                Kyu{i,2} = simwhiteXwhiteKernCompute(simwhiteKern1, whiteKern, x, x2);
            end
            Kuu = cell(2,1);
            Kuu{1} = rbfKernCompute(rbfKern, x2);
            Kuu{2} = noise*eye(length(x2));
            startVal = 1;
            endVal = 0;
            KuuMat = zeros(2*length(x2));
            for k = 1:nin,
                endVal = endVal + length(x2);
                KuuMat(startVal:endVal,startVal:endVal) = Kuu{k};
                startVal = endVal + 1;
            end
            K = [KuuMat cell2mat(Kyu)';cell2mat(Kyu) cell2mat(Kyy)];
            yu = real(gsamp(zeros(size(K,1),1), K, 1));
            u = yu(1:2*size(x2,1));
            y = yu(2*size(x2,1)+1:end);
            U = reshape(u,size(x2,1),nin);
            Y = reshape(y,size(x,1),nout);
            for k=1:nout,
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            ntrainx =100;
            maxl = length(x);
            X = cell(1,nout+nin);
            y = cell(1,nout+nin);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = x(pindx,:);
                y{k} = Y(pindx,k);
                yTest{k} = Y(indx(ntrainx+1:end),k);
            end
            XTest = x(indx(ntrainx+1:end),:)';
            % Append the latent functions
            y{nout+1} = U(:,1);
            y{nout+2} = U(:,2);
            X{nout+1} = x2;
            X{nout+2} = x2;
            save([baseDir 'simToyCombined.mat'], 'X', 'y', 'XTest', 'yTest')
        end
    case 'ggwhiteToyHighDimBatch'
        %         try
        %             load([baseDir 'ggwhiteToyHighDimBatch.mat']);
        %         catch
        kernName = 'ggwhite';
        inverseWidth = [10 5 5 100 30 200 400 20 1 40];
        sensitivity =  [1  2  2.2 1.5 1  3   4   2  1    10];
        noisePerOutput = 1e-2;
        nlf = 1;
        nout = 5;
        d = nlf + nout;
        inputDim = 3;
        % Sample the inputs
        X1 = gsamp(zeros(inputDim,1), eye(inputDim), 100);
        %            X1 = linspace(-1,1, 200)';
        kernType{1} = multigpKernComposer(kernName, d, nlf, 'ftc', 1);
        kernType{2} = multigpKernComposer('white',  d, nlf, 'ftc', 1);
        kern = kernCreate(X1,  {'cmpnd', kernType{:}});
        kern.comp{1}.comp{1}.variance = 1;
        for k = 1:nout,
            kern.comp{1}.comp{1+k}.precisionG = inverseWidth(k)*ones(inputDim,1);
            kern.comp{1}.comp{1+k}.variance = sensitivity(k);
            kern.comp{2}.comp{1+k}.variance = noisePerOutput;
        end
        K = kernCompute(kern, X1);
        yu = gsamp(zeros(size(K,1),1), K, 1);
        u = yu(1:size(X1,1));
        y = yu(size(X1,1)+1:end);
        Kout = K(1+size(X1,1):end,1+size(X1,1):end);
        [invKout, U, jitter] = pdinv(Kout);
        if any(jitter>1e-4)
            fprintf('Warning: Added jitter of %2.4f\n', jitter)
        end
        logDetKout = logdet(Kout, U);
        dim = size(y, 2);
        ll = -dim*log(2*pi) -logDetKout - y*invKout*y';
        ll = ll*0.5;
        U = reshape(u,size(X1,1),nlf);
        Y = reshape(y,size(X1,1),nout);
        X = cell(1, nout);
        y = cell(1, nout);
        for k =1:nout,
            X{k} = X1;
            y{k} = Y(:,k);
        end
        % Append the latent functions
        y{nout+1} = U(:,1);
        X{nout+1} = X1;
        XTest = [];
        yTest = [];
        save([baseDir 'ggwhiteToyHighDimBatch.mat'], 'X', 'y', 'XTest', 'yTest')
        %end

    case  'ggToyConvolution'
        try
            load([baseDir 'toyMultigp1DConvExample.mat']);
        catch
            randn('seed', 1e5);
            rand('seed', 1e5);
            nout = 5;
            nlf = 1;
            d = nout + nlf;
            precisionG = [10000 500 5000 100 1];
            variance = 5+5*rand(1,nout);
            precisionU = 500;
            sigma2Latent = 1;
            N = 500;
            N2 = 100;
            x = linspace(-1,1,N)';
            x2 = linspace(-1,1,N2)';
            XX = cell(1, d);
            XX{1} = x2;
            for i =1: nout;
                XX{nlf+i} = x;
            end
            kernType{1} = multigpKernComposer('gg', d, nlf, 'ftc', 1);
            kern = kernCreate(XX,  kernType{:});
            kern.comp{1}.precisionU = precisionU;
            kern.comp{1}.sigma2Latent = sigma2Latent;
            for k = 1:nout,
                kern.comp{1+k}.precisionU = precisionU;
                kern.comp{1+k}.sigma2Latent = sigma2Latent;
                kern.comp{1+k}.precisionG = precisionG(k);
                kern.comp{1+k}.sensitivity = variance(k);
            end
            K = kernCompute(kern, XX);
            yu = gsamp(zeros(size(K,1),1), K, 1);
            y = yu(size(XX{1},1)+1:end);
            Y = reshape(y,size(XX{2},1),nout);
            F = reshape(y,size(XX{2},1),nout);
            for k=1:nout,
                F(:,k) = Y(:,k);
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            X = cell(1,nout);
            y = cell(1,nout);
            f = cell(1,nout);
            for k =1:nout,
                X{k} = XX{k+nlf};
                y{k} = Y(:,k);
                f{k} = F(:,k);
            end
            XTest = [];
            yTest = f;
            save([baseDir 'toyMultigp1DConvExample.mat'], 'X', 'y', 'XTest', 'yTest');
        end

    case 'ggToyConvolutionTrainTest'
        try
            load([baseDir 'toyMultigp1DConvExampleTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        catch
            [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData('ggToyConvolution');
            nout = size(XTemp,2);
            ntrainx = 200;
            maxl = length(XTemp{1});
            X = cell(1,nout);
            y = cell(1,nout);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = XTemp{k}(pindx,:);
                y{k} = yTemp{k}(pindx,:);
                yTest{k} = yTemp{k}(indx(ntrainx+1:end),:);
            end
            XTest = XTemp{1}(indx(ntrainx+1:end),:)';
            yTest{end+1} = yTestTemp; %It has the groundtruth
            save([baseDir 'toyMultigp1DConvExampleTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        end

    case  'ggToyICM'
        try
            load([baseDir 'toyMultigp1DICMExample.mat']);
        catch
            randn('seed', 1e5);
            rand('seed', 1e5);
            nout = 5;
            nlf = 1;
            d = nout + nlf;
            rank = 5;
            Am = rand(nout,rank);
            precisionU = 50;
            sigma2Latent = 1;
            options.kern.nout = nout;
            options.kern.rankCorregMatrix = rank;
            N = 500;
            x = linspace(-1,1,N)';
            kernType{1} = multigpKernComposer('lmc', d, nlf, 'ftc', 1, options);
            kern = kernCreate(x,  kernType{:});
            kern.comp{1}.precisionU = precisionU;
            kern.comp{1}.sigma2Latent = sigma2Latent;
            kern.comp{2}.precisionU = precisionU;
            kern.comp{2}.sigma2Latent = sigma2Latent;
            kern.comp{2}.A = Am;
            kern.comp{2}.B = kern.comp{2}.A*kern.comp{2}.A';
            K = kernCompute(kern.comp{2}, x);
            y = gsamp(zeros(size(K,1),1), K, 1);
            Y = reshape(y,size(x,1),nout);
            F = reshape(y,size(x,1),nout);
            for k=1:nout,
                F(:,k) = Y(:,k);
                Y(:,k) = Y(:,k) + 0.1*sqrt(var(Y(:,k)))*randn(size(Y(:,k),1),1);
            end
            X = cell(1,nout);
            y = cell(1,nout);
            f = cell(1,nout);
            for k =1:nout,
                X{k} = x;
                y{k} = Y(:,k);
                f{k} = F(:,k);
            end
            XTest = [];
            yTest = f;
            figure
            for k=1:nout
                subplot(nout,1,k);
                plot(x,f{k})
            end
            save([baseDir 'toyMultigp1DICMExample.mat'], 'X', 'y', 'XTest', 'yTest');
        end

    case 'ggToyICMTrainTest'
        try
            load([baseDir 'toyMultigp1DICMTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        catch
            [XTemp, yTemp, XTestTemp, yTestTemp] = mapLoadData('ggToyICM');
            nout = size(XTemp,2);
            ntrainx = 200;
            maxl = length(XTemp{1});
            X = cell(1,nout);
            y = cell(1,nout);
            yTest = cell(1,nout);
            indx = randperm(maxl);
            pindx = sort(indx(1:ntrainx));
            for k =1:nout,
                X{k} = XTemp{k}(pindx,:);
                y{k} = yTemp{k}(pindx,:);
                yTest{k} = yTemp{k}(indx(ntrainx+1:end),:);
            end
            XTest = XTemp{1}(indx(ntrainx+1:end),:)';
            yTest{end+1} = yTestTemp; %It has the groundtruth
            save([baseDir 'toyMultigp1DICMTrainTest.mat'], 'X', 'y', 'XTest', 'yTest')
        end

    case 'sarcosMultiGP'
        dataTrain = load([baseDir 'sarcos_inv.mat']);
        dataTest = load([baseDir 'sarcos_inv_test.mat']);
        N = 2000;
        Ntest = 1000;
        trainD = dataTrain.sarcos_inv;
        testD = dataTest.sarcos_inv_test;
        sizeTrain = size(trainD);
        sizeTest = size(testD);
        Xmean = mean(trainD(:,1:21),1);
        Xstd = std(trainD(:,1:21),1);
        indexTraining = randperm(sizeTrain(1));
        indexTraining = indexTraining(1:N);
        indexTesting = randperm(sizeTest(1));
        indexTesting = indexTesting(1:Ntest);
        %   nTrainD = (trainD(:,1:21) - repmat(Xmean, sizeTrain(1),1 ))./ ...
        %           repmat(Xstd, sizeTrain(1),1);
        %   X = cell(1,7);
        %   y = cell(1,7);
        %   XTest = cell(1,7);
        %   yTest = cell(1,7);
        for i=1:7
            %       X{i} = (trainD(indexTraining,1:21) - repmat(Xmean, N,1 ))./ ...
            %           repmat(Xstd, N,1 );
            %       y{i} = trainD(indexTraining,21+i);
            %       XTest{i} = (testD(indexTesting,1:21) - repmat(Xmean, Ntest,1 ))./ ...
            %           repmat(Xstd, Ntest,1 );
            %       yTest{i} = testD(indexTesting,21+i);
            X{i} = trainD(indexTraining,1:21);
            %       y{i} = trainD(indexTraining,21+i);
            XTest{i} = testD(indexTesting,1:21);
            %       yTest{i} = testD(indexTesting,21+i);
        end
        X{1} = trainD(indexTraining,1:21);
        XTest{1} = testD(indexTesting,1:21);
        y{1} = trainD(indexTraining,22);
        yTest{1} = testD(indexTesting,22);

        save([baseDir 'sarcosData.mat'], 'X', 'y', 'XTest', 'yTest')

    case 'compilerData'
        data = load([baseDir 'data_compiler_org.mat']);
        X = cell(1,length(data.X));
        y = cell(1,length(data.X));
        %          for k = 1:length(data.X),
        %              X{k} = zscore(data.X{k});
        %              y{k} = zscore(data.Y{k});
        %          end
        X = data.X';
        y= data.Y';
        XTest = [];
        yTest = [];

    case 'robotWireless'
        [strength, time, x1, x2, storedMacs] = parseWirelessData([baseDir 'uw-floor.txt']);
        Xfull = [x1, x2];
        for i = 1:size(strength, 1)
            %[void, order] = sort(strength(i, :));
            %strength(i, order(1:end-8)) = NaN;
            strength(i, find(strength(i, :)<=-90)) = NaN;
        end
        y1 = strength(1:215, :);
        X1 = Xfull(1:215, :);

        yTest = cell(1, size(y1, 2));
        XTest = cell(1, size(y1, 2));
        y = cell(1, size(y1, 2));
        X = cell(1, size(y1, 2));

        yTest1 = strength(216:end, :);
        XTest1 = Xfull(216:end, :);
        indTrain  = zeros(1,1);
        contTrain = 0;
        indTest   = zeros(1,1);
        contTest  = 0;
        for i = 1:size(y1, 2)
            ind = find(~isnan(y1(:, i)));
            if isempty(ind) %|| length(ind)<=2
                contTrain = contTrain + 1;
                indTrain(contTrain) = i;
            end
            y{i} = y1(ind, i);
            X{i} = X1(ind, :);
            ind = find(~isnan(yTest1(:, i)));
            if isempty(ind) %|| length(ind)<=2
                contTest = contTest + 1;
                indTest(contTest) = i;
            end
            yTest{i} = yTest1(ind, i);
            XTest{i} = XTest1(ind, :);
        end
        X([indTrain indTest]) = [];
        y([indTrain indTest]) = [];
        XTest([indTrain indTest]) = [];
        yTest([indTrain indTest]) = [];

    case 'schoolData1'
        % This characterization of the School Data is due to Bakker and
        % Keskes
        try
            load([baseDir 'schoolData1.mat']);
        catch
            nout = 139;
            file = [baseDir 'ILEA567.DAT'];
            fid = fopen(file,'r');
            cont = 0;
            nStudents = 15362;
            % Student-dependent features
            yearExam = zeros(nStudents, 3);          %3 features (dummy variables)
            gender = zeros(nStudents, 1);              %1 features (dummy variables)
            Vrband = zeros(nStudents, 1);              % 1 features (dummy variables)
            ethnicGroup = zeros(nStudents, 11);   % 11 features (dummy variables)
            % School-dependent features
            perEligibleStudents = zeros(nStudents, 1); % 1 feature
            VR1band = zeros(nStudents, 1); % 1 feature
            schoolGender = zeros(nStudents, 1); % 1 features
            schoolDenomination = zeros(nStudents, 3); % 3 features (dummy variables)
            % Task
            task = zeros(nStudents, 1); %
            % Exam score
            examScore = zeros(nStudents, 1);
            for k =1:nStudents
                cont = cont + 1;
                rawData = fgetl(fid);
                % STUDENT-Dependent Feature
                % Year of exam
                index = str2double(rawData(1));
                yearExam(cont, index) = 1;
                % Gender
                value = str2double(rawData(11));
                gender(cont, 1) = value;
                % VR band of student
                index = str2double(rawData(12));
                Vrband(cont, 1) = index;
                % Ethnic group
                index = str2double(rawData(13:14));
                ethnicGroup(cont, index) = 1;
                % SCHOOL-Dependent Features
                % Percent. Students eligible for free school meals
                value = str2double(rawData(7:8));
                perEligibleStudents(cont, 1) = value;
                % Percent. Students in school in VR band 1
                value = str2double(rawData(9:10));
                VR1band(cont, 1) = value;
                % School gender
                index = str2double(rawData(15));
                if index == 1;
                    schoolGender(cont, 1) = 0;
                else
                    schoolGender(cont, 1) = 1;
                end
                % School denomination
                index = str2double(rawData(16));
                schoolDenomination(cont, index) = 1;
                % TASK index
                value = str2double(rawData(2:4));
                task(cont, 1) = value;
                % EXAM scores == OUTPUTS
                value = str2double(rawData(5:6));
                examScore(cont, 1) = value;
                indicator = feof(fid);
                if indicator
                    break;
                end
            end
            % Organize the tasks and the inputs per task
            features = [yearExam gender Vrband ethnicGroup ...
                perEligibleStudents VR1band schoolGender schoolDenomination];
            fclose(fid);
            features =  features(:, 1:16);
            % Normalization as in Bakker and Heskes et al
            features = zscore(features);
            examScore = zscore(examScore);
            X = cell(1,nout);
            y = cell(1,nout);
            XTest = cell(1,nout);
            yTest = cell(1,nout);
            for j=1:nout,
                %                 XT = features(task == j, :);
                %                 XT(:, 1:16) = zscore(XT(:, 1:16));
                %                 X{j}  = XT;
                %                 y{j}= zscore(examScore(task == j, :));
                X{j}  = features(task == j, :);
                y{j}= examScore(task == j, :);
            end
            save([baseDir 'schoolData1.mat'],'X','y','XTest', 'yTest');
        end

    case 'schoolData2'
        % This characterization of the School data is due to Bonilla et al
        try
            load([baseDir 'schoolData2.mat']);
        catch
            nout = 139;
            file = [baseDir 'ILEA567.DAT'];
            fid = fopen(file,'r');
            cont = 0;
            nStudents = 15362;
            % Student-dependent features
            yearExam = zeros(nStudents, 3); % 3 features (dummy variables)
            gender = zeros(nStudents, 2); % 1 features (dummy variables)
            Vrband = zeros(nStudents, 4); % 4 features (dummy variables) (Michelli)
            ethnicGroup = zeros(nStudents, 11); % 11 features (dummy variables)
            % School-dependent features
            perEligibleStudents = zeros(nStudents, 1); % 1 feature
            VR1band = zeros(nStudents, 1); % 1 feature
            schoolGender = zeros(nStudents, 3); % 3 features (dummy)
            schoolDenomination = zeros(nStudents, 3); % 3 features (dummy variables)
            % Task
            task = zeros(nStudents, 1); %
            % Exam score
            examScore = zeros(nStudents, 1);
            for k =1:nStudents
                cont = cont + 1;
                rawData = fgetl(fid);
                % STUDENT-Dependent Feature
                % Year of exam
                index = str2double(rawData(1));
                yearExam(cont, index) = 1;
                % Gender
                value = str2double(rawData(11));
                gender(cont, 1+ value) = 1;
                % VR band of student
                index = str2double(rawData(12));
                Vrband(cont, index + 1) = 1;
                % Ethnic group
                index = str2double(rawData(13:14));
                ethnicGroup(cont, index) = 1;
                % SCHOOL-Dependent Features
                % Percent. Students eligible for free school meals
                value = str2double(rawData(7:8));
                perEligibleStudents(cont, 1) = value;
                % Percent. Students in school in VR band 1
                value = str2double(rawData(9:10));
                VR1band(cont, 1) = value;
                % School gender
                index = str2double(rawData(15));
                schoolGender(cont, index) = 1;
                % School denomination
                index = str2double(rawData(16));
                schoolDenomination(cont, index) = 1;
                % TASK index
                value = str2double(rawData(2:4));
                task(cont, 1) = value;
                % EXAM scores == OUTPUTS
                value = str2double(rawData(5:6));
                examScore(cont, 1) = value;
                indicator = feof(fid);
                if indicator
                    break;
                end
            end
            % Organize the tasks and the inputs per task
            features = [yearExam gender Vrband ethnicGroup ...
                perEligibleStudents VR1band schoolGender schoolDenomination];
            fclose(fid);
            features = features(:,1:20); % Only student dependent features
            features = zscore(features); % Normaliza the features
            XTemp = cell(nout,1);
            yTemp = cell(nout,1);
            for j=1:nout,
                XTemp{j} = features(task == j, :);
                yTemp{j} = examScore(task == j, :);
            end
            X = cell(1,nout);
            y = cell(1,nout);
            XTest = cell(1,nout);
            yTest = cell(1,nout);
            % This bit finds unique features as in Bonilla et
            % al paper. Just to test what is best.
            cont = 0;
            q = size(features,2);
            for j=1:nout,
                [uniqueX, I, J] = unique(XTemp{j},'rows');
                [sorted, indexJ] = sort(J);
                for k=1: size(uniqueX,1),
                    indexToAvg = find(sorted == k);
                    X{j}(k,1:q) = XTemp{j}(indexJ(indexToAvg(1)), :);
                    % We append the number of vectors per output point.
                    X{j}(k,q+1) = length(indexToAvg);
                    y{j}(k,:) = mean(yTemp{j}(indexJ(indexToAvg), 1));
                end
                cont = cont + size(uniqueX, 1);
            end
            save([baseDir 'schoolData2.mat'],'X','y','XTest', 'yTest');
        end

    case 'schoolData3'
        % This characterization of the School data is due to Michelli et al
        try
            load([baseDir 'schoolData3.mat']);
        catch
            nout = 139;
            file = [baseDir 'ILEA567.DAT'];
            fid = fopen(file,'r');
            cont = 0;
            nStudents = 15362;
            % Student-dependent features
            yearExam = zeros(nStudents, 3); % 3 features (dummy variables)
            gender = zeros(nStudents, 1); % 1 features (dummy variables)
            Vrband = zeros(nStudents, 4); % 4 features (dummy variables) (Michelli)
            ethnicGroup = zeros(nStudents, 11); % 11 features (dummy variables)
            % School-dependent features
            perEligibleStudents = zeros(nStudents, 1); % 1 feature
            VR1band = zeros(nStudents, 1); % 1 feature
            schoolGender = zeros(nStudents, 3); % 3 features (dummy)
            schoolDenomination = zeros(nStudents, 3); % 3 features (dummy variables)
            % Task
            task = zeros(nStudents, 1); %
            % Exam score
            examScore = zeros(nStudents, 1);
            for k =1:nStudents
                cont = cont + 1;
                rawData = fgetl(fid);
                % STUDENT-Dependent Feature
                % Year of exam
                index = str2double(rawData(1));
                yearExam(cont, index) = 1;
                % Gender
                value = str2double(rawData(11));
                gender(cont, 1) = value;
                % VR band of student
                index = str2double(rawData(12));
                %Vrband(cont, index + 1) = 1;
                Vrband(cont, 1) = index;
                % Ethnic group
                index = str2double(rawData(13:14));
                ethnicGroup(cont, index) = 1;
                % SCHOOL-Dependent Features
                % Percent. Students eligible for free school meals
                value = str2double(rawData(7:8));
                perEligibleStudents(cont, 1) = value;
                % Percent. Students in school in VR band 1
                value = str2double(rawData(9:10));
                VR1band(cont, 1) = value;
                % School gender
                schoolGender(cont, index) = 1;
                % School denomination
                index = str2double(rawData(16));
                schoolDenomination(cont, index) = 1;
                % TASK index
                value = str2double(rawData(2:4));
                task(cont, 1) = value;
                % EXAM scores == OUTPUTS
                value = str2double(rawData(5:6));
                examScore(cont, 1) = value;
                indicator = feof(fid);
                if indicator
                    break;
                end
            end
            % Organize the tasks and the inputs per task
            features = [yearExam gender Vrband ethnicGroup ...
                perEligibleStudents VR1band schoolGender schoolDenomination];
            fclose(fid);
            % Normalization of the features.
            features = zscore(features);
            X = cell(1,nout);
            y = cell(1,nout);
            XTest = cell(1,nout);
            yTest = cell(1,nout);
            for j=1:nout,
                X{j} = features(task == j, :);
                y{j} = examScore(task == j, :);
            end
            save([baseDir 'schoolData3.mat'],'X','y','XTest', 'yTest');
        end

    case 'juraDataCd'
        try
            load([baseDir 'juraDataCd.mat']);
        catch
            fidP = fopen([baseDir 'prediction.dat'],'r');
            if fidP ==-1
                error('The file prediction.dat does not exist in this directory');
            end
            fidV = fopen([baseDir 'validation.dat'],'r');
            if fidV ==-1
                error('The file validation.dat does not exist in this directory');
            end
            xLoc = [1 2];
            primaryLoc = 5; % Column location of Cd in file prediction.dat and file validation.dat
            secondaryLoc = [9 11]; % Column location of Ni and Zn
            fseek(fidP, 61,-1);
            A = textscan(fidP, '%f%f%f%f%f%f%f%f%f%f%f');
            predValues = cell2mat(A);
            fseek(fidV, 61,-1);
            A = textscan(fidV, '%f%f%f%f%f%f%f%f%f%f%f');
            valValues = cell2mat(A);
            X = cell(1, 1+length(secondaryLoc));
            y = cell(1, 1+length(secondaryLoc));
            XTest =  cell(1, 1+length(secondaryLoc));
            yTest =  cell(1, 1+length(secondaryLoc));
            X{1} = predValues(:, xLoc);
            y{1} = predValues(:,primaryLoc);
            for i=1:length(secondaryLoc);
                X{1+i} = predValues(:, xLoc);
                y{1+i} = predValues(:, secondaryLoc(i));
            end
            % Append validation values to training data
            for i=1:length(secondaryLoc);
                X{1+i} = [X{1+i}; valValues(:, xLoc)];
                y{1+i} = [y{1+i}; valValues(:, secondaryLoc(i))];
            end
            % Form the testting sets
            XTest{1} = valValues(:, xLoc);
            yTest{1} = valValues(:, primaryLoc);
            for i=1:length(secondaryLoc);
                XTest{1+i} = valValues(:, xLoc);
                yTest{1+i} = valValues(:, secondaryLoc(i));
            end
            fclose(fidP);
            fclose(fidV);
            save([baseDir 'juraDataCd.mat'], 'X', 'y', 'XTest', 'yTest');

        end
        %/~
    case 'juraDataCo'
        data = load([baseDir 'data_jura_Co']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;
        %~/
    case 'juraDataCu'
        try
            load([baseDir 'juraDataCu.mat']);
        catch
            fidP = fopen([baseDir 'prediction.dat'],'r');
            if fidP ==-1
                error('The file prediction.dat does not exist in this directory');
            end
            fidV = fopen([baseDir 'validation.dat'],'r');
            if fidV ==-1
                error('The file validation.dat does not exist in this directory');
            end
            xLoc = [1 2];
            primaryLoc = 8; % Column location of Cu in file prediction.dat and file validation.dat
            secondaryLoc = [10 9 11]; % Column location of Ni and Zn
            fseek(fidP, 61,-1);
            A = textscan(fidP, '%f%f%f%f%f%f%f%f%f%f%f');
            predValues = cell2mat(A);
            fseek(fidV, 61,-1);
            A = textscan(fidV, '%f%f%f%f%f%f%f%f%f%f%f');
            valValues = cell2mat(A);
            X = cell(1, 1+length(secondaryLoc));
            y = cell(1, 1+length(secondaryLoc));
            XTest =  cell(1, 1+length(secondaryLoc));
            yTest =  cell(1, 1+length(secondaryLoc));
            X{1} = predValues(:, xLoc);
            y{1} = predValues(:,primaryLoc);
            for i=1:length(secondaryLoc);
                X{1+i} = predValues(:, xLoc);
                y{1+i} = predValues(:, secondaryLoc(i));
            end
            % Append validation values to training data
            for i=1:length(secondaryLoc);
                X{1+i} = [X{1+i}; valValues(:, xLoc)];
                y{1+i} = [y{1+i}; valValues(:, secondaryLoc(i))];
            end
            % Form the testting sets
            XTest{1} = valValues(:, xLoc);
            yTest{1} = valValues(:, primaryLoc);
            for i=1:length(secondaryLoc);
                XTest{1+i} = valValues(:, xLoc);
                yTest{1+i} = valValues(:, secondaryLoc(i));
            end
            fclose(fidP);
            fclose(fidV);
            save([baseDir 'juraDataCu.mat'], 'X', 'y', 'XTest', 'yTest');
        end
        %/~
    case 'juraDataPb'
        data = load([baseDir 'data_jura_Pb']);
        X = data.Xtrain;
        y = data.Ytrain;
        XTest = data.Xtest;
        yTest = data.Ytest;

    case 'juraDataCdIsotopic'
        try
            load([baseDir 'juraDataCdIsotopic.mat']);
        catch
            fidP = fopen([baseDir 'prediction.dat'],'r');
            if fidP ==-1
                error('The file prediction.dat does not exist in this directory');
            end
            fidV = fopen([baseDir 'validation.dat'],'r');
            if fidV ==-1
                error('The file validation.dat does not exist in this directory');
            end
            xLoc = [1 2];
            primaryLoc = 5; % Column location of Cd in file prediction.dat and file validation.dat
            secondaryLoc = [9 11]; % Column location of Ni and Zn
            fseek(fidP, 61,-1);
            A = textscan(fidP, '%f%f%f%f%f%f%f%f%f%f%f');
            predValues = cell2mat(A);
            fseek(fidV, 61,-1);
            A = textscan(fidV, '%f%f%f%f%f%f%f%f%f%f%f');
            valValues = cell2mat(A);
            X = cell(1, 1+length(secondaryLoc));
            y = cell(1, 1+length(secondaryLoc));
            XTest =  cell(1, 1+length(secondaryLoc));
            yTest =  cell(1, 1+length(secondaryLoc));
            X{1} = predValues(:, xLoc);
            y{1} = predValues(:,primaryLoc);
            for i=1:length(secondaryLoc);
                X{1+i} = predValues(:, xLoc);
                y{1+i} = predValues(:, secondaryLoc(i));
            end
            %     % Append validation values to training data
            %     for i=1:length(secondaryLoc);
            %       X{1+i} = [X{1+i}; valValues(:, xLoc)];
            %       y{1+i} = [y{1+i}; valValues(:, secondaryLoc(i))];
            %     end
            % Form the testting sets
            XTest{1} = valValues(:, xLoc);
            yTest{1} = valValues(:, primaryLoc);
            for i=1:length(secondaryLoc);
                XTest{1+i} = valValues(:, xLoc);
                yTest{1+i} = valValues(:, secondaryLoc(i));
            end
            fclose(fidP);
            fclose(fidV);
            save([baseDir 'juraDataCdIsotopic.mat'], 'X', 'y', 'XTest', 'yTest');

        end
    
        
    case 'demp53_5genes'
        data = load([baseDir 'dataBarencoOption_0_Genes_5']);
        X = data.data{1}.Xtrain;
        y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
        XTest = X;
        %yTest = y;
        % We use the variance
        yTest = [cell2mat(data.data{1}.Yvar); cell2mat(data.data{2}.Yvar); cell2mat(data.data{3}.Yvar)];

    case 'demp53_50genes'
        data = load([baseDir 'dataBarencoOption_0_Genes_50']);
        X = data.data{1}.Xtrain;
        y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
        XTest = X;
        % yTest = y;
        % We use the variance
        yTest = [cell2mat(data.data{1}.Yvar); cell2mat(data.data{2}.Yvar); cell2mat(data.data{3}.Yvar)];

    case 'demp53_20genes'
        data = load([baseDir 'dataBarencoOption_0_Genes_20']);
        X = data.data{1}.Xtrain;
        y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
        XTest = X;
        %yTest = y;
        % We use the variance
        yTest = [cell2mat(data.data{1}.Yvar); cell2mat(data.data{2}.Yvar); cell2mat(data.data{3}.Yvar)];
        %~/

    case 'demp53_54genes'
        data = load([baseDir 'dataBarencoOption_0_Genes_54']);
        X = data.data{1}.Xtrain;
        y = [cell2mat(data.data{1}.Ytrain); cell2mat(data.data{2}.Ytrain); cell2mat(data.data{3}.Ytrain)];
        XTest = X;
        %yTest = y;
        % We use the variance
        yTest = [cell2mat(data.data{1}.Yvar); cell2mat(data.data{2}.Yvar); cell2mat(data.data{3}.Yvar)];

    case 'demDrosOnlyTwiRed'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 50;
        maxGenes = size(Genes, 2);
        indexes = randperm(maxGenes);
        X = timePoints;
        y = Genes(:, indexes(1:nGenesRed));
        XTest = X;
        yTest = GenesVar(:, indexes(1:nGenesRed));

    case 'demDrosOnlyTwiRed100'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 100;
        maxGenes = size(Genes, 2);
        indexes = randperm(maxGenes);
        X = timePoints;
        y = Genes(:, indexes(1:nGenesRed));
        XTest = X;
        yTest = GenesVar(:, indexes(1:nGenesRed));

    case 'demDrosOnlyTwi'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 1500;
        maxGenes = size(Genes, 2);
        indexes = randperm(maxGenes);
        X = timePoints;
        y = Genes(:, indexes(1:nGenesRed));
        XTest = X;
        yTest = GenesVar(:, indexes(1:nGenesRed));

    case 'demDros92'

        load([baseDir 'dros92']);
        %   nGenesRed = 1500;
        %   maxGenes = size(Genes, 2);
        %   indexes = randperm(maxGenes);
        X = timePoints;
        y = Genes;
        XTest = X;
        yTest = GenesVar;


    case 'demDrosOnlyTwiFiltered1'
        load([baseDir 'demDrosOnlyTwiFiltered1']);

    case 'demDrosOnlyTwiFiltered2'
        load([baseDir 'demDrosOnlyTwiFiltered2']);

    case 'demDrosARed'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 50;
        maxGenes = size(Genes, 2);
        indexes = randperm(maxGenes);
        indexes = randperm(maxGenes);
        indexes = randperm(maxGenes);
        yTemp = Genes(:, indexes(1:nGenesRed));
        indexR1 = 1:12;
        indexR3 = 25:36;
        y = cell(nGenesRed,1);
        yTest = cell(nGenesRed, 1);
        X = cell(nGenesRed,1);
        XTest = cell(nGenesRed, 1);
        for i=1:nGenesRed
            y{i} = yTemp(indexR1, i);
            yTest{i} = yTemp(indexR3, i);
            X{i} = timePoints;
            XTest{i} = timePoints;
        end

    case 'demDrosBRed'
        load([baseDir 'drosOnlyTwist']);
        load([baseDir 'intersecDrosRBF']);
        yTemp = Genes(:, indexesFilter);
        nGenesRed = 50;
        indexes = randperm(length(intersec13));
        subIntersec13 = indexes(1:nGenesRed);
        indexR1 = 1:12;
        indexR3 = 25:36;
        y = cell(length(subIntersec13),1);
        yTest = cell(length(subIntersec13), 1);
        X = cell(length(subIntersec13),1);
        XTest = cell(length(subIntersec13), 1);
        for i=1:length(subIntersec13)
            y{i} = yTemp(indexR1, subIntersec13(i));
            yTest{i} = yTemp(indexR3, subIntersec13(i));
            X{i} = timePoints;
            XTest{i} = timePoints;
        end

    case 'demDrosCRed'
        load([baseDir 'drosOnlyTwist']);
        load([baseDir 'intersecDrosSIM']);
        yTemp = Genes(:, indexesFilter);
        nGenesRed = 50;
        indexes = randperm(length(intersec13));
        subIntersec13 = indexes(1:nGenesRed);
        indexR1 = 1:12;
        indexR3 = 25:36;
        y = cell(length(subIntersec13),1);
        yTest = cell(length(subIntersec13), 1);
        X = cell(length(subIntersec13),1);
        XTest = cell(length(subIntersec13), 1);
        for i=1:length(subIntersec13)
            y{i} = yTemp(indexR1, subIntersec13(i));
            yTest{i} = yTemp(indexR3, subIntersec13(i));
            X{i} = timePoints;
            XTest{i} = timePoints;
        end

    case 'demDrosA'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 1000;
        yTemp = Genes(:, 1:nGenesRed);
        indexR1 = 1:12;
        indexR3 = 25:36;
        y = cell(nGenesRed,1);
        yTest = cell(nGenesRed, 1);
        X = cell(nGenesRed,1);
        XTest = cell(nGenesRed, 1);
        for i=1:nGenesRed
            y{i} = yTemp(indexR1, i);
            yTest{i} = yTemp(indexR3, i);
            X{i} = timePoints;
            XTest{i} = timePoints;
        end

    case 'demDrosAComplete'

        load([baseDir 'drosOnlyTwist']);
        nGenesRed = 1500;
        yTemp = Genes(:, 1:nGenesRed);
        indexR1 = 1:12;
        indexR3 = 25:36;
        y = cell(nGenesRed,1);
        yTest = cell(nGenesRed, 1);
        X = cell(nGenesRed,1);
        XTest = cell(nGenesRed, 1);
        for i=1:nGenesRed
            y{i} = yTemp(indexR1, i);
            yTest{i} = yTemp(indexR3, i);
            X{i} = timePoints;
            XTest{i} = timePoints;
        end


    case 'yeastSpellman'
        load([baseDir 'yeastSpellman']);
        yTest = X; % Connectivity matrix
        X = (1:size(data,2))';
        y = data;
        %
        XTest = [];

    case 'yeastSpellmanRed'
        load([baseDir 'yeastSpellmanRed']);
        yTest = X; % Connectivity matrix
        X = (1:size(data,2))';
        y = data;
        %
        XTest = [];

    case 'sensorsTemperature'

        data = load([baseDir 'temperatureSensor.mat']);
        X = data.XTrain;
        y = data.yTrain;
        XTest = data.XTest;
        yTest = data.yTest;

    case 'toySdlfm1'
        try
            load([baseDir 'toySdlfm1.mat']);
        catch
            Xbase = linspace(0, 20, 200)';
            numOut = 1;
            nSamples = 10;
            X = cell(1, numOut);
            for i=1:numOut
                X{i} = Xbase;
            end
            options.nIntervals = 3;
            options.nlfPerInt = 1;
            options.approx = 'ftc';
            options.includeVel = false;
            options.includeAccel = false;
            options.kernType = 'sdlfm';
            options.kern.includeVel = options.includeVel;
            options.kern.includeAccel = options.includeAccel;
            options.kern.nIntervals = options.nIntervals;
            options.kern.isNormalised = true;
            kernType = cell(1);
            F = cell(1, nSamples);
            Y = cell(1, nSamples);
            for q=1:nSamples
                options.kern.switchingTimes = [-1  6+0.1*randn 7+0.1*randn];
                options.kern.nlfPerInt = options.nlfPerInt;
                kernType{1} = sdlfmgpKernComposer(options.kernType, numOut, 'ftc', options);
                kern = kernCreate(X, {'cmpnd', kernType{:}});
                kern = sdmultiKernParamInit(kern, options);
                params = kernExtractParam(kern);
                ind = paramNameRegularExpressionLookup(kern, '.* inverse width .*', true);
                params(ind) = log([1e-3 1e-3 1e-3]);
                ind = paramNameRegularExpressionLookup(kern, '.* mass', true);
                params(ind) = log(0.1);
                ind = paramNameRegularExpressionLookup(kern, '.* damper', true);
                params(ind) = log(0.4);
                ind = paramNameRegularExpressionLookup(kern, '.* spring', true);
                params(ind) = log(5);
                ind = paramNameRegularExpressionLookup(kern, '.* sensitivity .*', true);
                params(ind) =  [10 20 10];
                kern = kernExpandParam(kern, params);
                %%% Pass the parameters to the sdrbf kernel
                % We need first to create the SDRBF kernel and pass the parameters
                % to it.
                sdrbfKern = kernCreate(Xbase, {'parametric', options.kern, 'sdrbf'});
                indxIW = paramNameRegularExpressionLookup(kern, '.* inverse width .*',1);
                indxSP = paramNameRegularExpressionLookup(kern, '.* switching point .*',1);
                indxIW1 = paramNameRegularExpressionLookup(sdrbfKern, 'inverse width .*',1);
                indxSP1 = paramNameRegularExpressionLookup(sdrbfKern, 'switching point .*',1);
                params = kernExtractParam(kern);
                params1 = kernExtractParam(sdrbfKern);
                params1(indxIW1) = params(indxIW);
                params1(indxSP1) = params(indxSP);
                sdrbfKern = kernExpandParam(sdrbfKern, params1);
                Kfu = zeros(numOut*length(Xbase), options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:numOut
                    endVal = endVal + length(X{i});
                    KTemp = sdlfmXsdrbfKernCompute(kern.comp{1}.comp{i}, sdrbfKern, X{i}, Xbase);
                    Kfu(startVal:endVal, :) = cell2mat(KTemp);
                    startVal = endVal + 1;
                end
                KuuTemp = sdrbfKernCompute(sdrbfKern, Xbase);
                KuuM = zeros(options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:options.nlfPerInt
                    endVal = endVal + length(Xbase);
                    KuuM(startVal:endVal, startVal:endVal) = KuuTemp{i};
                    startVal = endVal + 1;
                end
                Kff = kernCompute(kern, X);
                K = [KuuM  Kfu';Kfu Kff];
                f = real(gsamp(zeros(size(K,1),1), K, 1))';
                f = f(length(Xbase)+1:end,:);
                F{q} = f;
                Y{q} = f + 0.1*sqrt(var(f))*randn(length(Xbase),1);
            end
            XTest = X;
            y = Y;
            yTest = F;
            save([baseDir 'toySdlfm1.mat'], 'X', 'y', 'XTest', 'yTest');
        end
    case 'toySdlfm2'
        try
            load([baseDir 'toySdlfm2.mat']);
        catch
            Xbase = linspace(0, 20, 500)';
            numOut = 2;
            nSamples = 10;
            X = cell(1, numOut);
            for i=1:numOut
                X{i} = Xbase;
            end
            options.nIntervals = 3;
            options.nlfPerInt = 1;
            options.approx = 'ftc';
            options.includeVel = false;
            options.includeAccel = false;
            options.kernType = 'sdlfm';
            options.kern.includeVel = options.includeVel;
            options.kern.includeAccel = options.includeAccel;
            options.kern.nIntervals = options.nIntervals;
            options.kern.isNormalised = true;
            mass = 0.1;
            damper = [0.4 1];
            spring = [2 3];
            inverseWidth = 1e-3*ones(1,options.nIntervals);
            sensitivity = [10 1 10 5 -10 1];
            meanSwitching = [6 7];
            numberTraining = 200;
            kernType = cell(1);
            F = cell(1, nSamples);
            Y = cell(1, nSamples);
            Ytest = cell(1, nSamples);
            Ytrain = cell(1, nSamples);
            Xtest = cell(1, nSamples);
            Xtrain = cell(1, nSamples);
            for q=1:nSamples
                options.kern.switchingTimes = [-1  meanSwitching(1) meanSwitching(2)];
                options.kern.nlfPerInt = options.nlfPerInt;
                kernType{1} = sdlfmgpKernComposer(options.kernType, numOut, 'ftc', options);
                kern = kernCreate(X, {'cmpnd', kernType{:}});
                kern = sdmultiKernParamInit(kern, options);
                params = kernExtractParam(kern);
                ind = paramNameRegularExpressionLookup(kern, '.* mass', true);
                params(ind) = log(mass);
                ind = paramNameRegularExpressionLookup(kern, '.* damper', true);
                params(ind) = log(damper);
                ind = paramNameRegularExpressionLookup(kern, '.* spring', true);
                params(ind) = log(spring);
                ind = paramNameRegularExpressionLookup(kern, '.* sensitivity .*', true);
                params(ind) =  sensitivity;
                kern = kernExpandParam(kern, params);
                for i=1:numOut
                    kern.comp{1}.comp{i}.inverseWidth = inverseWidth;
                end
                %%% Pass the parameters to the sdrbf kernel
                % We need first to create the SDRBF kernel and pass the parameters
                % to it.
                sdrbfKern = kernCreate(Xbase, {'parametric', options.kern, 'sdrbf'});
                sdrbfKern.inverseWidth = inverseWidth;
                Kfu = zeros(numOut*length(Xbase), options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:numOut
                    endVal = endVal + length(X{i});
                    KTemp = sdlfmXsdrbfKernCompute(kern.comp{1}.comp{i}, sdrbfKern, X{i}, Xbase);
                    Kfu(startVal:endVal, :) = cell2mat(KTemp);
                    startVal = endVal + 1;
                end
                KuuTemp = sdrbfKernCompute(sdrbfKern, Xbase);
                KuuM = zeros(options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:options.nlfPerInt
                    endVal = endVal + length(Xbase);
                    KuuM(startVal:endVal, startVal:endVal) = KuuTemp{i};
                    startVal = endVal + 1;
                end
                Kff = kernCompute(kern, X);
                K = [KuuM  Kfu';Kfu Kff];
                f = real(gsamp(zeros(size(K,1),1), K, 1))';
                f = f(length(Xbase)+1:end,:);
                startVal = 1;
                endVal = 0;
                indexShuffled = randperm(length(Xbase));
                indexTrain = sort(indexShuffled(1:numberTraining));
                indexTest = sort(indexShuffled(numberTraining+1:end));
                for i=1:numOut
                    endVal = endVal + length(Xbase);
                    F{q}{i} = f(startVal:endVal);
                    Y{q}{i} = F{q}{i} + 0.05*sqrt(var(F{q}{i}))*randn(length(Xbase),1);
                    Ytrain{q}{i} = Y{q}{i}(indexTrain);
                    Ytest{q}{i} = Y{q}{i}(indexTest);
                    Xtrain{q}{i} = Xbase(indexTrain);
                    Xtest{q}{i} = Xbase(indexTest);
                    startVal = endVal + 1;
                end
            end
            y = Ytrain;
            X = Xtrain;
            XTest = Xtest;
            XTest{end+1} = Xbase';
            yTest = Ytest;
            yTest{end+1} = F;
            save([baseDir 'toySdlfm2.mat'], 'X', 'y', 'XTest', 'yTest');
        end
    case 'toySdlfm3'
        try
            load([baseDir 'toySdlfm3.mat']);
        catch
            Xbase = linspace(0, 20, 500)';
            numOut = 3;
            nSamples = 10;
            X = cell(1, numOut);
            for i=1:numOut
                X{i} = Xbase;
            end
            options.nIntervals = 2;
            options.nlfPerInt = 1;
            options.approx = 'ftc';
            options.includeVel = false;
            options.includeAccel = false;
            options.kernType = 'sdlfm';
            options.kern.includeVel = options.includeVel;
            options.kern.includeAccel = options.includeAccel;
            options.kern.nIntervals = options.nIntervals;
            options.kern.isNormalised = true;
            mass = 0.1;
            damper = [0.4 1 1];
            spring = [2 3 0.5];
            inverseWidth = [1e-3 1];
            sensitivity = [1 5 5 1 1 1];
            meanSwitching = 10;
            numberTraining = 200;
            kernType = cell(1);
            F = cell(1, nSamples);
            Y = cell(1, nSamples);
            Ytest = cell(1, nSamples);
            Ytrain = cell(1, nSamples);
            Xtest = cell(1, nSamples);
            Xtrain = cell(1, nSamples);
            for q=1:nSamples
                options.kern.switchingTimes = [-2  meanSwitching(1)];
                options.kern.nlfPerInt = options.nlfPerInt;
                kernType{1} = sdlfmgpKernComposer(options.kernType, numOut, 'ftc', options);
                kern = kernCreate(X, {'cmpnd', kernType{:}});
                kern = sdmultiKernParamInit(kern, options);
                params = kernExtractParam(kern);
                ind = paramNameRegularExpressionLookup(kern, '.* mass', true);
                params(ind) = log(mass);
                ind = paramNameRegularExpressionLookup(kern, '.* damper', true);
                params(ind) = log(damper);
                ind = paramNameRegularExpressionLookup(kern, '.* spring', true);
                params(ind) = log(spring);
                ind = paramNameRegularExpressionLookup(kern, '.* sensitivity .*', true);
                params(ind) =  sensitivity;
                kern = kernExpandParam(kern, params);
                for i=1:numOut
                    kern.comp{1}.comp{i}.inverseWidth = inverseWidth;
                end
                %%% Pass the parameters to the sdrbf kernel
                % We need first to create the SDRBF kernel and pass the parameters
                % to it.
                sdrbfKern = kernCreate(Xbase, {'parametric', options.kern, 'sdrbf'});
                sdrbfKern.inverseWidth = inverseWidth;
                Kfu = zeros(numOut*length(Xbase), options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:numOut
                    endVal = endVal + length(X{i});
                    KTemp = sdlfmXsdrbfKernCompute(kern.comp{1}.comp{i}, sdrbfKern, X{i}, Xbase);
                    Kfu(startVal:endVal, :) = cell2mat(KTemp);
                    startVal = endVal + 1;
                end
                KuuTemp = sdrbfKernCompute(sdrbfKern, Xbase);
                KuuM = zeros(options.nlfPerInt*length(Xbase));
                startVal = 1;
                endVal = 0;
                for i=1:options.nlfPerInt
                    endVal = endVal + length(Xbase);
                    KuuM(startVal:endVal, startVal:endVal) = KuuTemp{i};
                    startVal = endVal + 1;
                end
                Kff = kernCompute(kern, X);
                K = [KuuM  Kfu';Kfu Kff];
                f = real(gsamp(zeros(size(K,1),1), K, 1))';
                f = f(length(Xbase)+1:end,:);
                startVal = 1;
                endVal = 0;
                indexShuffled = randperm(length(Xbase));
                indexTrain = sort(indexShuffled(1:numberTraining));
                indexTest = sort(indexShuffled(numberTraining+1:end));
                for i=1:numOut
                    endVal = endVal + length(Xbase);
                    F{q}{i} = f(startVal:endVal);
                    Y{q}{i} = F{q}{i} + 0.1*sqrt(var(F{q}{i}))*randn(length(Xbase),1);
                    Ytrain{q}{i} = Y{q}{i}(indexTrain);
                    Ytest{q}{i} = Y{q}{i}(indexTest);
                    Xtrain{q}{i} = Xbase(indexTrain);
                    Xtest{q}{i} = Xbase(indexTest);
                    startVal = endVal + 1;
                end
            end
            y = Ytrain;
            X = Xtrain;
            XTest = Xtest;
            XTest{end+1} = Xbase';
            yTest = Ytest;
            yTest{end+1} = F;
            save([baseDir 'toySdlfm3.mat'], 'X', 'y', 'XTest', 'yTest');
        end

    case 'demDrosMel'
        try
            load([baseDir 'drosMelPLOS.mat']);
        catch
            fidraw = fopen([baseDir 'journal.pcbi.0020051.sd007.txt'],'r');
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'w');
            indicator = feof(fidraw);
            while ~indicator
                currentLine = fgetl(fidraw);
                index = regexp(currentLine, '\d');
                if ~isempty(index)
                    if length(index) > 5
                        fprintf(fiddata, '%s\n', currentLine);
                    end
                end
                indicator = feof(fidraw);
            end
            fclose(fidraw);
            fclose(fiddata);
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'r');
            A = textscan(fiddata, '%f%f%f%f%f%f%f');
            fclose(fiddata);
            timePoints = [0;11;24;30;37;43;49;55;62;68];
            spacePoints = 0.5*(0.35:0.01:0.92)';
            y = cell(1, 7);
            X = cell(1, 7);
            yTest = cell(1, 7);
            XTest = cell(1, 7);
            % Create input space
            XM = zeros(size(A{1},1),2);
            startOne = 1;
            endOne = 0;
            for i=1:length(timePoints)
                endOne = endOne + 58;
                XM(startOne:endOne, 1) = timePoints(i);
                XM(startOne:endOne, 2) = spacePoints;
                startOne = endOne + 1;
            end
            for i=1:7
                y{i} = A{i};
                X{i} = XM;
            end
            yTest = y;
            XTest = X;
            save([baseDir 'drosMelPLOS.mat'], 'X', 'y', 'XTest', 'yTest');
        end

    case 'demDrosMelInf'
        try
            load([baseDir 'drosMelPLOSInf.mat']);
        catch
            fidraw = fopen([baseDir 'journal.pcbi.0020051.sd007.txt'],'r');
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'w');
            indicator = feof(fidraw);
            while ~indicator
                currentLine = fgetl(fidraw);
                index = regexp(currentLine, '\d');
                if ~isempty(index)
                    if length(index) > 5
                        fprintf(fiddata, '%s\n', currentLine);
                    end
                end
                indicator = feof(fidraw);
            end
            fclose(fidraw);
            fclose(fiddata);
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'r');
            A = textscan(fiddata, '%f%f%f%f%f%f%f');
            fclose(fiddata);
            timePoints = [0;11;24;30;37;43;49;55;62;68];
            spacePoints = 0.5*(0.35:0.01:0.92)';
            y = cell(1, 7);
            X = cell(1, 7);
            yTest = cell(1, 7);
            XTest = cell(1, 7);
            % Create input space
            XM = zeros(size(A{1},1),2);            
            XM(1:length(timePoints), 1) = timePoints;
            XM(1:length(spacePoints), 2) = spacePoints;
            XM(length(timePoints)+1:580,1) = Inf;
            XM(length(spacePoints)+1:580,2) = Inf;
            for i=1:7
                y{i} = A{i};
                X{i} = XM;
            end
            yTest = y;
            XTest = X;
            save([baseDir 'drosMelPLOSInf.mat'], 'X', 'y', 'XTest', 'yTest');
        end
        
    case 'demDrosMelRed'
        try
            load([baseDir 'drosMelPLOSRed.mat']);
        catch
            fidraw = fopen([baseDir 'journal.pcbi.0020051.sd007.txt'],'r');
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'w');
            indicator = feof(fidraw);
            while ~indicator
                currentLine = fgetl(fidraw);
                index = regexp(currentLine, '\d');
                if ~isempty(index)
                    if length(index) > 5
                        fprintf(fiddata, '%s\n', currentLine);
                    end
                end
                indicator = feof(fidraw);
            end
            fclose(fidraw);
            fclose(fiddata);
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'r');
            A = textscan(fiddata, '%f%f%f%f%f%f%f');
            fclose(fiddata);
            timePoints = [0;11;24;30;37;43;49;55;62;68];
            spacePoints = 0.5*(0.35:0.01:0.92)';
            lengthInit = size(A{1},1);
            samplingRate = 4;
            redA = cell(length(timePoints),1);
            for i=1:7
                vecA = mat2cell(A{i}, length(spacePoints)*ones(length(timePoints),1), 1);
                for j=1:length(timePoints)
                    redA{j} = vecA{j}(1:samplingRate:length(spacePoints));                    
                end
                A{i} = cell2mat(redA);
            end
            spacePoints = spacePoints(1:samplingRate:length(spacePoints));
            y = cell(1, 7);
            X = cell(1, 7);
            yTest = cell(1, 7);
            XTest = cell(1, 7);
            % Create input space
            XM = zeros(size(A{1},1),2);
            startOne = 1;
            endOne = 0;
            for i=1:length(timePoints)
                endOne = endOne + length(spacePoints);
                XM(startOne:endOne, 1) = timePoints(i);
                XM(startOne:endOne, 2) = spacePoints;
                startOne = endOne + 1;
            end
            for i=1:7
                y{i} = A{i};
                X{i} = XM;
            end
            yTest = y;
            XTest = X;
            save([baseDir 'drosMelPLOSRed.mat'], 'X', 'y', 'XTest', 'yTest');
        end
        
    case 'demDrosMelRedInf'
        try
            load([baseDir 'drosMelPLOSRedInf.mat']);
        catch
            fidraw = fopen([baseDir 'journal.pcbi.0020051.sd007.txt'],'r');
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'w');
            indicator = feof(fidraw);
            while ~indicator
                currentLine = fgetl(fidraw);
                index = regexp(currentLine, '\d');
                if ~isempty(index)
                    if length(index) > 5
                        fprintf(fiddata, '%s\n', currentLine);
                    end
                end
                indicator = feof(fidraw);
            end
            fclose(fidraw);
            fclose(fiddata);
            fiddata = fopen([baseDir 'drosMelPLOS.txt'],'r');
            A = textscan(fiddata, '%f%f%f%f%f%f%f');
            fclose(fiddata);
            timePoints = [0;11;24;30;37;43;49;55;62;68];
            spacePoints = 0.5*(0.35:0.01:0.92)';
            lengthInit = size(A{1},1);
            samplingRate = 4;
            redA = cell(length(timePoints),1);
            for i=1:7
                vecA = mat2cell(A{i}, length(spacePoints)*ones(length(timePoints),1), 1);
                for j=1:length(timePoints)
                    redA{j} = vecA{j}(1:samplingRate:length(spacePoints));
                end
                A{i} = cell2mat(redA);
            end
            spacePoints = spacePoints(1:samplingRate:length(spacePoints));
            y = cell(1, 7);
            X = cell(1, 7);
            % Create input space
            XM = zeros(size(A{1},1),2);
            XM(1:length(timePoints), 1) = timePoints;
            XM(1:length(spacePoints), 2) = spacePoints;
            XM(length(timePoints)+1:end,1) = Inf;
            XM(length(spacePoints)+1:end,2) = Inf;
            for i=1:7
                y{i} = A{i};
                X{i} = XM;
            end
            yTest = y;
            XTest = X;
            save([baseDir 'drosMelPLOSRedInf.mat'], 'X', 'y', 'XTest', 'yTest');
        end
    otherwise
        error('Unknown data set requested.')

end

if semiSup % Test if data is for semi-supervised learning.
    indUnlabelled = find(rand(size(y, 1), 1)>labProb);
    y(indUnlabelled, :) = NaN;
end
