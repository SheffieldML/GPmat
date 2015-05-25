% CLASSIFICATIONDEMO2 Try the IVM for classification.

noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
display = 2;
dVal = 200;
% Sample a regression data-set.
generateClassificationData;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)
for i = 1:4
  if display > 1
    clf
    pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
    pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
    set(pointsNeg, 'erasemode', 'xor')
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, 100);
end

model = ivmOptimiseIVM(model, display);


% model.kernelType = 'ARD';
% prior = 0;
% display = 2;

% % Sample a classification data-set.
% generateClassificationData;
% model.X = X;
% model.y = y;

% numData = size(model.X, 1);
% switch model.kernelType
%  case 'regular'
%   model.lntheta = log(10*[ones(1, 5)]);
%  case 'ARD'
%   model.lntheta = log([10 10 1 1 1 ones(1, size(model.X, 2))]);
% end


% % Number of inclusions
% model.d = 1000;


% % selection criteria
% model.selectionCriteria = 'entropy';
% model.noiseType = 'probit';
% nClass1 = sum(model.y==1);
% nClass2 = sum(model.y==-1);
% model.bias = invCummGaussian(nClass1/(nClass2+nClass1));
% %model.noiseVariance = 0.001;
clf
plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
hold on
plot(X(find(y==1), 1), X(find(y==1), 2), 'ro')
%model = ivmOptimiseIVM(model, display);

