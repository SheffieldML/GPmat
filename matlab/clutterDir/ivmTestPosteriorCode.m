% IVMTESTPOSTERIORCODE Test the posterior code.

randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,3)-[zeros(100, 1) 6*ones(100, 2)]; randn(100,3)+[zeros(100, 1) 6*ones(100, 2)]; randn(100, 3)];
y = [ones(200, 1); -ones(100, 1)];

noiseModel = 'gaussian';
selectionCriterion = 'entropy';
kernelType = 'regular';
prior = 0;
display = 1;
dVal = 10;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)
for i = 1:0
  if display > 1
    cl
    pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
    set(pointsNeg, 'erasemode', 'xor')
    hold on
    pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
    set(pointsNeg, 'erasemode', 'xor')
  end
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, 100);
end
if display > 1
  clf
  pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
  set(pointsNeg, 'erasemode', 'xor')
  hold on
  pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
  set(pointsNeg, 'erasemode', 'xor')
end
model = ivmOptimiseIVM(model, display);

origX = randn(1, 3);
epsilon = 1e-6;

for i = 1:3
  X = origX;
  X(i) = origX(i)+epsilon;
  [dmu(i), dvarsigma(i)] = ivmPosteriorMeanVar(model, X);
  fhandle = str2func([model.noise.type 'LogLikelihood']);
  L(i) = fhandle(X, 1, model);
  X(i) = origX(i)-epsilon;
  [dmum(i), dvarsigmam(i)] = ivmPosteriorMeanVar(model, X);
  Lm(i) = fhandle(X, 1, model);

end
gmu = 0.5*(dmu - dmum)'/epsilon;
gvarsigma = 0.5*(dvarsigma - dvarsigmam)'/epsilon;
gL = 0.5*(L - Lm)/epsilon;
[gmu2, gvarsigma2] = ivmPosteriorGradMeanVar(model, X);
fhandle = str2func([model.noise.type 'GradX']);
gL2 =  fhandle(X, 1, model);