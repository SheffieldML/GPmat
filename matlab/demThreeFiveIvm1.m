% DEMTHREEFIVEIVM1 Try the IVM & NCNM on 3 vs 5.

% IVM


rand('seed', 1e5);
randn('seed', 1e5);

dataSetName = 'threeFive';
experimentNo = 1;

[X, y, XTest, yTest] = mapLoadData(dataSetName);

% IVM Set up

optionsIvm = ivmOptions;
dVal = 200;
optionsIvm.display = 0;
optionsIvm.kernIters = 400;
optionsIvm.extIters = 15;
optionsIvm.kern = {'rbf', 'white'};

% SVM Set up
optionsSvm = svmlopt;
optionsSvm.Kernel = 2;


% Use this prior for NCNM Process variances.
prior.type = 'gamma';
prior = priorParamInit(prior);
prior.a = 1;
prior.b = 1;
num = 6;
probLabelled = [0.2 0.1 0.05 0.025 0.0125];

for i = 1:10
  for j = 1:length(probLabelled)
    rand('seed', i*1e5);
    randn('seed', i*1e5);

    prob = probLabelled(j);
    unlabelled = find(rand(size(y))>prob);
    
    % Create labelled only data-sets.
    ylab = y;
    Xlab = X;
    ylab(unlabelled) = [];
    Xlab(unlabelled, :) = [];
    yunlab = y;
    % Run IVM on labelled data.
    optionsIvm.noise = 'probit';
    optionsIvm.numActive = min([dVal length(ylab)]);

    % Initialise the model.
    model = ivmCreate(size(Xlab, 1), size(ylab, 2), Xlab, ylab, optionsIvm);

    model.noise.bias = 0;
    % Optimise the IVM but don't change noise params.
    model = ivmOptimise(model, optionsIvm)
    % Select data-points in an IVM with those kernel parameters.
    model = ivmOptimiseIvm(model, optionsIvm.display);
    % Display the final model.
    ivmDisplay(model);
    [mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
    areaIvm(i, j) = rocCurve(mu, yTest);
    fprintf('IVM %d, %d: Area under ROC curve: %2.4f\n', i, j, areaIvm(i, j));
    % Store the results.
    [kernIvm(i, j), noiseIvm(i, j), ivmInfoIvm(i, j)] = ivmDeconstruct(model);
    
    % Run NCNM Ivm on data.
    optionsIvm.noise = 'ncnm';
    optionsIvm.numActive = dVal;
    yunlab(unlabelled) = NaN;
    % Intitalise NCNM
    model = ivmCreate(size(X, 1), size(yunlab, 2), X, yunlab, optionsIvm);
    model.noise.bias = 0;
    model.noise.width = 1;
    prior.index = 2;
    model.kern.comp{1}.priors(1) = prior;
    prior.index = 1;
    model.kern.comp{2}.priors(1) = prior;
    % Optimise the NCNM but don't change noise params.
    model = ivmOptimise(model, optionsIvm);
    % Select data-points in an NCNM with those kernel parameters.
    model = ivmOptimiseIvm(model, optionsIvm.display);
    % Display the final model.
    ivmDisplay(model);
    [mu, varSigma] = ivmPosteriorMeanVar(model, XTest);
    areaNCNM(i, j) = rocCurve(mu, yTest);
    fprintf('NCNM %d, %d: Area under ROC curve: %2.4f\n', i, j, areaNCNM(i, j));
    % Store the results.
    [kernNCNM(i, j), noiseNCNM(i, j), ivmInfoNCNM(i, j)] = ivmDeconstruct(model);
    
    % Run SVM on labelled data.
    optionsSvm.KernelParam = kernIvm(i, j).comp{1}.inverseWidth/2;
    model = svml([dataSetName num2str(i) num2str(j)], optionsSvm);
    model = svmltrain(model, Xlab, ylab);
    [alpha, Xactive] = svmlread(model.fname);
    % Fix bug in svmlread which doesn't recognise size of Xactive
    missingCols = size(XTest, 2) - size(Xactive, 2);
    Xactive = [Xactive zeros(size(Xactive, 1), missingCols)];
    missingRows = size(alpha, 1) - size(Xactive, 1);
    Xactive = [Xactive; zeros(missingRows, size(Xactive, 2))];
    % Make output predictions.
    testOut = XTest*Xactive'*alpha;
    areaSVM(i, j) = rocCurve(testOut, yTest);
    fprintf('SVM %d, %d, Area under ROC curve: %2.4f\n', i, j, areaSVM(i, j));
    
    % Run Transductive SVM on data.
    yunlab(unlabelled) = 0;
    optionsSvm.KernelParam = max([kernNCNM(i, j).comp{1}.inverseWidth kernIvm(i, j).comp{1}.inverseWidth])/2;
    model = svml([dataSetName num2str(i) num2str(j)], optionsSvm);
    model = svmltrain(model, X, yunlab);
    [alpha, Xactive] = svmlread(model.fname);
    % Fix bug in svmlread which doesn't recognise size of Xactive
    missingCols = size(XTest, 2) - size(Xactive, 2);
    Xactive = [Xactive zeros(size(Xactive, 1), missingCols)];
    missingRows = size(alpha, 1) - size(Xactive, 1);
    Xactive = [Xactive; zeros(missingRows, size(Xactive, 2))];
    % Make output predictions.
    testOut = XTest*Xactive'*alpha;
    areaTSVM(i, j) = rocCurve(testOut, yTest);
    fprintf('T-SVM %d, %d, Area under ROC curve: %2.4f\n', i, j, areaTSVM(i, j));
    capName = dataSetName;
    capName(1) = upper(capName(1));
    save(['dem' capName 'Ivm1.mat'], ...
         'kernIvm', 'noiseIvm', ...
         'ivmInfoIvm', 'areaIvm', ...
         'kernNCNM', 'noiseNCNM', ...
         'ivmInfoNCNM', 'areaNCNM', ...
         'probLabelled', 'areaSVM', 'areaTSVM');
  end
end


