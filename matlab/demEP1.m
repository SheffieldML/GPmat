% DEMEP1 Demonstrate Expectation propagation.


randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
y = [ones(200, 1); -ones(100, 1)];

noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
display = 2;
dVal = 30;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)
if display > 1
  clf
  pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
  set(pointsNeg, 'erasemode', 'xor')
  hold on
  pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
  set(pointsNeg, 'erasemode', 'xor')
end
model = ivmOptimiseIVM(model, display);
b = model.beta;
counter = 0;

% while 1
%   counter = counter + 1;
%   index = randperm(length(model.I));
%   for j = index
%     i = model.I(j);
%     model = ivmEPUpdate(model, i);
%   end
%   change = max(abs(model.beta - b));
%   fprintf('EP iteration %d, diff %2.4f\n', counter, full(change));
%   if  change < 1e-6
%     break
%   end
%   if ~rem(counter, 10)
%     fprintf('Paused\n')
%     pause
%   end
% end
% if display > 1
%   clf
%   pointsNeg = plot(X(find(y==-1), 1), X(find(y==-1), 2), 'bx');
%   set(pointsNeg, 'erasemode', 'xor')
%   hold on
%   pointsPos = plot(X(find(y==1), 1), X(find(y==1), 2), 'ro');
%   set(pointsNeg, 'erasemode', 'xor')
% end





