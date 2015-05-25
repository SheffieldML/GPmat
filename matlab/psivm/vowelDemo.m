% CLASSIFICATIONDEMO Try the IVM for classification.

model.kern.type = 'ARD';
prior = 0;
display = 1;

% Sample a classification data-set.
generateVowelData;
model.X = X;
model.y = vowelY(:, 2);

numData = size(model.X, 1);
switch model.kern.type
 case 'regular'
  model.kern.lntheta = log(10*[ones(1, 5)]);
 case 'ARD'
  model.kern.lntheta = log([10 10 1 1 1 ones(1, size(model.X, 2))]);
end


% Number of inclusions
model.d = 200;


% selection criteria
model.selectionCriteria = 'entropy';
model.noise.type = 'probit';
nClass1 = sum(model.y==1);
nClass2 = sum(model.y==-1);
model.noise.bias = invCummGaussian(nClass1/(nClass2+nClass1));
%model.noiseVariance = 0.001;
for i = 1:4

  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior);
  disp(exp(model.kern.lntheta))
end
model = ivmOptimiseIVM(model, display);

