% DEMREGRESSIONSIMPLE Demonstrate Gaussian processes for regression.
% FORMAT
% DESC runs a simple one-D Gaussian process displaying errorbars.
%
% SEEALSO : demRegression, gpCreate, demInterpolation
%
% COPYRIGHT: Andreas C. Damianou, 2011, 2012
% COPYRIGHT: Neil Lawrence 2006

% SHEFFIELDML

% Uncomment to fix seeds so that every run will reproduce the same results
%randn('seed', 1e7)
%rand('seed', 1e7)

close all


% The hyperparameters of the GP covariance function (kernel) can be
% optimised with some gradient-based optimiser, or can just be set by hand
% to fixed values
optimiseHyperparameters = false;
noiseLevel = 0.1;
% Assume noise, i.e. in the generative model, y's are generated from x's
% according to : y = f(x) + e, where f ~GP and e is Gaussian noise.
% If noiseLevel=0, then in the predictions, the uncertainty around observed
% (training) points will be zero.
noiseVar = noiseLevel*noiseLevel;

%----- Create data set (x, yTrue)
x = linspace(-1, 1, 50)'; % input indices in the x-axis
trueKern = kernCreate(x, 'rbf'); % Generate samples from this kernel
trueKern.inverseWidth = 10;
K = kernCompute(trueKern, x) + eye(size(x, 1))*noiseVar;
% Sample some true function values.
yTrue = gsamp(zeros(size(x))', K, 1)';
plot(x, yTrue, 'x-'), title('True function'), figure
%%
%----- Fitting a GP ---------
%... by assuming that only a subset of the (x,yTrue)
% pairs are observed
permInd = randperm(size(x,1));
percTrInputs = 0.5;
indTrain = permInd(1:floor(percTrInputs*size(x,1)))';
indTest = setdiff(1:size(x,1), indTrain)';
%indTrain = [1 3 5 8 9]'; % Experiment with different sub-vectors of [1:9]
%indTrain = [1:9]';


% (xTrain, yTrain) is the subset of (x,yTrue) that is presented to the
% model
yTrain = yTrue(indTrain);
xTrain = x(indTrain);

% Test inputs:
% We want to get values from the fitted GP in the interval -2:2 with large
% density (so that it looks like a continuous line)
xTest = linspace(-2, 2, 200)';



% Standard GP prediction:

if ~optimiseHyperparameters
    % No optimisation for kernel hyperparameters and no sparsity
    % approximation. Here we just use the standard equations of GP
    % prediction.
    
    % We wish to fit our data using an RBF kernel (doesn't necessarily have to
    % be the same kernel that generated the data, e.g. experiment with
    % 'matern32' etc. In theory, we shouldn't even know how the data were
    % generated. See the KERN toolbox for all kernels and their hyperparameters).
    kern = kernCreate(x, 'rbf');
    % Change inverse variance (1/(lengthScale^2)))
    % Set manually the hyperparameter for RBF. In a real-world application, the
    % hyperparameters can be learned via optimisation (see above).
    kern.inverseWidth = 2;
    
    % Evaluate the cov. matrix which contains also the cross-covariance
    % between xTest and xTrain
    Kx = kernCompute(kern, xTest, xTrain);
    % The covariance matrix evaluated only on the training points
    Ktrain = kernCompute(kern, xTrain, xTrain);
    
    % The means of the prediction
    yPred = Kx*pdinv(Ktrain + eye(size(Ktrain))*noiseVar)*yTrain;
    % Variance of the prediction
    yVar = kernDiagCompute(kern, xTest) - sum(Kx*pdinv(Ktrain+ eye(size(Ktrain))*noiseVar).*Kx, 2);
else
    % Optimise hyperparameters. This option allows to also experiment with
    % the sparse GP framework, by giving to gpOptions the appropriate
    % argument, depending on the desired approximation scheme.
    
    optionsGP = gpOptions('ftc');
    
    % Uncomment the following to obtain a sparse GP approximation (faster
    % but less accurate)
    %optionsGP = gpOptions('fitc');
    %optionsGP.numActive = floor(size(xTrain,1)/2);
    
    % The kernel used is a sum of an rbf, a bias and a white term. Bias and
    % white help numerical stability during optimisation (this will also
    % result in bigger error-bars in the end). Also try 'lin' or 'matern32'
    % or even 'rbfperiodic' for the first component of the kernel.
    optionsGP.kern = {'matern32','bias','white'};
    %optionsGP.scale2var1 = true;     % Scale outputs to variance 1.
    modelGP = gpCreate(size(xTrain,2), size(yTrain,2), xTrain, yTrain, optionsGP);
    modelGP = gpOptimise(modelGP, 1, 1000);
    % The model is now optimised. The following function just does standard
    % GP prediction based on the optimised model and the test inputs.
    [yPred, yVar] = gpPosteriorMeanVar(modelGP, xTest);
end

ySd = sqrt(yVar);



%% -------- Plots---------
% These plots produce: a) A continuous line which is the GP prediction on
% the xTest vector which is very dense so that the line looks continuous
% b) training input - output pairs denoted with 'x'
% c) test input - output pairs denoted with 'o'. These are the REAL test
% values. The GP prediction (continuous line) should nicely interpolate
% between them (because we KNOW that we observe training data under noise)
% and not be expected to pass exactly from the true test values. In fact,
% the predicted continuous line should not EXACTLY pass from all training
% values either, because in that case it would have overfitted (if the
% model is optimised).


% The gray area represents the uncertainty (due to the lack of training
% points nearby). The gray area is the errorbars, i.e. +- 2*Standard Dev.
fillColor =[0.7 0.7 0.7];
fill([xTest; xTest(end:-1:1)], ...
    [yPred; yPred(end:-1:1)] ...
    + 2*[ySd; -ySd], ...
    fillColor,'EdgeColor',fillColor)
hold on;
% Predictions (plot it as a continuous line to look like a function)
h = plot(xTest, yPred, 'k-');
hold on;
% Training points
p = plot(xTrain, yTrain, 'x', 'Color', 'r');

% The points that we considered to be unobserved
%xTestReal = setdiff(x,xTrain);
%yTestReal = setdiff(yTrue,yTrain);
xTestReal = x(indTest);
yTestReal = yTrue(indTest);
p2 = plot(xTestReal, yTestReal,'o');

% If the model is in sparse approximation mode, also plot the inducing
% points learned
if optimiseHyperparameters & ~strcmp(modelGP.approx, 'ftc')
    hold on; plot(modelGP.X_u, zeros(size(modelGP.X_u,1)), 'm+')
    legend('Uncertainty due to noise)', 'Predicted function values','training points', 'test points','inducing points')
else
    legend('Uncertainty due to noise)', 'Predicted function values','training points', 'test points')
end
title('GP regression')
xlabel('X');
ylabel('f(X)+noise');

%---- Quantify the error:
% Compute the mean squared error between the predicted values and the real
% ones, which correspond to the x-indices that were left out from the
% training set.
mserror = mean(abs(yTestReal(:) - yPred(setdiff(1:length(x), indTrain))));
fprintf('Mean Squared Error for the points that were considered to be unobserved (test points): %f \n', mserror);
