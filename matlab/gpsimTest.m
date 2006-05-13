% GPSIMTEST Test the single input motif code.

numSamp = 50;
t = linspace(0, 2, numSamp)';
interSpace = t(2)-t(1);
X_0_1 = 1;
C1 = 5;
D1 = 5;
delta1 = 0;
X_0_2 = 2;
C2 = 0.5;
D2 = 0.5;
delta2 = delta1;
sigma = 1;
rbfKern = kernCreate(t, 'rbf');
rbfKern.inverseWidth = 2/(sigma*sigma);
rbfKern.variance = 1;
simKern1 = kernCreate(t, 'sim');
simKern1.inverseWidth = 2/(sigma*sigma);
simKern1.variance = 1;
simKern2 = simKern1;
simKern1.initVal = X_0_1;
simKern2.initVal = X_0_2;
simKern1.variance = C1*C1;
simKern2.variance = C2*C2;
simKern1.decay = D1;
simKern2.decay = D2;
simKern1.delay = delta1;
dimKern2.delay = delta2;

% K_11 = kernCompute(rbfKern, t);
% K_22 = kernCompute(simKern1, t);
% K_33 = kernCompute(simKern2, t);
% K_21 = simXrbfKernCompute(simKern1, rbfKern, t);
% K_31 = simXrbfKernCompute(simKern2, rbfKern, t);
% K_12 = K_21';
% K_13 = K_31';
% K_23 = simXsimKernCompute(simKern1, simKern2, t);
% K_32 = K_23';

% K2 = [K_11 K_12 K_13; K_21 K_22 K_23; K_31 K_32 K_33];

kern = kernCreate(t, {'multi', 'rbf', 'sim', 'sim'});
kern.comp{1} = rbfKern;
kern.comp{2} = simKern1;
kern.comp{3} = simKern2;

K = kernCompute(kern, t);

x1 = -1;
while any(x1<0) | any(x2<0) | any(f<0)
  y = gsamp(zeros(1, size(K, 1)), K, 1);
  x1 = y(numSamp+1:2*numSamp)';
  x2 = y(2*numSamp+1:end)';
  f = y(1:numSamp)';
  x1 = x1 + X_0_1*exp(-D1*t);
  x2 = x2 + X_0_2*exp(-D2*t);
end
figure(1), clf, plot(f), hold on, plot(x1, 'c'), plot(x2, 'r')

f1Guess = (diff(x1)/interSpace+D1*(x1(1:end-1)))/C1;
f2Guess = (diff(x2)/interSpace+D2*(x2(1:end-1)))/C2;
figure(2), clf, plot(f), hold on, plot(f1Guess, 'c'), plot(f2Guess, 'r')