numPoints = 100;
numBumps = 5;
part{1} = linspace(0, 20, 20)'; % Middle finger
part{2} = linspace(20, 45, 20)'; % Index finger
part{3} = linspace(45, 75, 20)'; % Flap between finger and thumb
part{4} = linspace(75, 115, 20)'; % thumb
part{5} = linspace(115, 320, 30)'; % Palm
part{6} = linspace(320, 340, 20)'; % little finger
part{7} = linspace(340, 360, 20)'; % ring finger
theta = [];
for i = 1:length(part)
  theta = [theta; part{i}(1:end-1)];
end
theta = deg2rad(theta);
kern = kernCreate(theta, 'gibbsperiodic');

options = rbfperiodicOptions(numBumps);
kern.lengthScaleFunc = modelCreate('rbfperiodic', kern.inputDimension, 1, options);

sigma = [5 5 10 5 5];
centers = [10 32 90 330 350];

kern.lengthScaleFunc.sigma2 = deg2rad(sigma).*deg2rad(sigma);
kern.lengthScaleFunc.thetaBar = deg2rad(centers);
kern.lengthScaleFunc.bias = 1.5;
kern.lengthScaleFunc.weights = -1.5*ones(numBumps, 1);
kern.lengthScaleFunc.weights(3) = -1;
K = kernCompute(kern, theta);

for i = 1:16
  subplot(4, 4, i)
  y = gsamp(zeros(size(theta)), K, 2)';
  plot(y(:, 1), y(:, 2), 'rx-');
end
