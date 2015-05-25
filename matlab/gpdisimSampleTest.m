% GPDISIMSAMPLETEST Test the single input motif code.
% FORMAT
% DESC tests the GPDISIM kernels by sampling an f and x for a couple
% of genes. Double checks that everything is working by solving the
% differential equation numerically for each gene and comparing to
% f.

% GPDISIM

randn('seed', 1e5);
rand('seed', 1e5);

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
simKern = kernCreate(t, 'sim');
simKern.inverseWidth = 2/(sigma*sigma);
simKern.variance = 1;
simKern.decay = 0.1;
disimKern1 = kernCreate(t, 'disim');
disimKern1.inverseWidth = 2/(sigma*sigma);
disimKern1.variance = 1;
disimKern2 = disimKern1;
disimKern1.initVal = X_0_1;
disimKern2.initVal = X_0_2;
disimKern1.variance = C1*C1;
disimKern2.variance = C2*C2;
disimKern1.decay = D1;
disimKern2.decay = D2;
disimKern1.delay = delta1;
dimKern2.delay = delta2;

% K_11 = kernCompute(rbfKern, t);
% K_22 = kernCompute(disimKern1, t);
% K_33 = kernCompute(disimKern2, t);
% K_21 = disimXrbfKernCompute(disimKern1, rbfKern, t);
% K_31 = disimXrbfKernCompute(disimKern2, rbfKern, t);
% K_12 = K_21';
% K_13 = K_31';
% K_23 = disimXdisimKernCompute(disimKern1, disimKern2, t);
% K_32 = K_23';

% K2 = [K_11 K_12 K_13; K_21 K_22 K_23; K_31 K_32 K_33];

kern = kernCreate(t, {'multi', 'rbf', 'sim', 'disim', 'disim'});
kern.comp{1} = rbfKern;
kern.comp{2} = simKern;
kern.comp{3} = disimKern1;
kern.comp{4} = disimKern2;
counter = 0;
K = kernCompute(kern, t);
colordef white
counter = counter + 1;
figure(counter)
imagesc(K, [-1.1 1.1]);
handle = [];
set(gca, 'fontname', 'times')
handle = [handle; text('Interpreter', 'latex', 'string', '$$y(t)$$', 'position', ...
            [25 210])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$f(t)$$', 'position', ...
            [75 210])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_1(t)$$', 'position', ...
            [125 210])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_2(t)$$', 'position', ...
            [175 210])];
set(handle, 'horizontalalignment', 'right')
handle = [handle; text('Interpreter', 'latex', 'string', '$$y(t)$$', 'position', ...
            [-10 25])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$f(t)$$', 'position', ...
            [-10 75])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_1(t)$$', 'position', ...
            [-10 125])];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_2(t)$$', 'position', ...
            [-10 175])];
set(handle, 'horizontalalignment', 'center')

set(handle, 'fontsize', 34)
set(gca, 'fontsize', 34)
axis off
colorbar
print('-depsc', ['../tex/diagrams/gpdisimTestKernelImage'...
                 '.eps'])

% figure
% for i = 1:3
%   x1 = -1;
%   while any(x1<0) | any(x2<0) | any(f<0) | any(y<0)
%     s = gsamp(zeros(1, size(K, 1)), K+eye(size(K))*1e-3, 1);
%     f = s(numSamp+1:2*numSamp)';
%     x1 = s(2*numSamp+1:3*numSamp)';
%     x2 = s(3*numSamp+1:end)';
%     y = s(1:numSamp)';
%     x1 = x1 + X_0_1*exp(-D1*t);
%     x2 = x2 + X_0_2*exp(-D2*t);
%   end
  
%   counter = counter + 1;
%   figure(counter), clf, 
%   lhand = plot(f), hold on, lhand= [lhand plot(y, 'g')]; lhand = [lhand plot(x1, 'c')];, lhand = ...
%           [lhand plot(x2, 'r')];
%   set(lhand, 'linewidth', 4)
%   set(gca, 'fontsize', 18)
%   ylim = get(gca, 'ylim');
%   set(gca, 'ylim', ylim);
%   print('-depsc', ['../tex/diagrams/gpdisimTestSamples' num2str(counter) '.eps'])
  
%   f1Guess = (diff(x1)/interSpace+D1*(x1(1:end-1)))/C1;
%   f2Guess = (diff(x2)/interSpace+D2*(x2(1:end-1)))/C2;
%   counter = counter + 1;
%   figure(counter), clf
%   lhand = plot(f), hold on, lhand = [lhand plot(f1Guess, 'c')]; 
%   lhand = [lhand plot(f2Guess, 'r')];
%   set(lhand, 'linewidth', 4)
%   set(gca, 'fontsize', 18)
%   set(gca, 'ylim', ylim);
%   print('-depsc', ['../tex/diagrams/gpdisimTestSamples' num2str(counter) '.eps'])
% end