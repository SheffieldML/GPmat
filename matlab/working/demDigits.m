% DIGITSDEMO Try the IVM on some digits data.

% IVM

model.kernelType = 'ARD';
prior = 0;
display = 1;
load ../data/usps_train
fourIndex = find(ALL_T==4);
nineIndex = find(ALL_T==9);
model.selectionCriterion = 'entropy';
model.X = [ALL_DATA(fourIndex, :); ALL_DATA(nineIndex, :)];
model.y = [ones(length(fourIndex), 1); -ones(length(nineIndex), 1)];

numData = size(model.X, 1);
switch model.kernelType
 case 'regular'
  model.lntheta = [ones(1, 5)];
 case 'ARD'
  model.lntheta = [ones(1, 5) ones(1, size(model.X, 2))];
end


% Number of inclusions
model.d = 200;

nClass1 = sum(model.y==1);
nClass2 = sum(model.y== -1);

% Set the bias parameter
model.bias = invCummGaussian(nClass1/(nClass2+nClass1));

% selection criteria
model.selectionCriteria = 'entropy';
model.noiseType = 'probit';
for i = 1:4
 model = ivmOptimiseIVM(model, display);
 model = ivmOptimiseKernel(model, prior);
 disp(model.lntheta)
end
model.d = 200
model = ivmOptimiseIVM(model, display);

load c:\datasets\usps\matlab\16x16\usps_test
fourIndex = find(ALL_T==4);
nineIndex = find(ALL_T==9);

xTest = [ALL_DATA(fourIndex, :); ALL_DATA(nineIndex, :)];
yTest = [ones(length(fourIndex), 1); -ones(length(nineIndex), 1)];
yPred = ivmfwd(xTest, model);

testError = 1-sum(yPred==yTest)/size(yTest, 1);