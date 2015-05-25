% ICMLTOYSINE A small demo of the multi-task IVM.

% PSIVM

%clear all
%close all

seed = 1e4;
prior = 0;
display = 0;
innerIters = 100; % Number of scg iterations
outerIters = 4;

kernelType = 'rbf';
noiseType = 'gaussian';
selectionCriterion = 'entropy';



d = 15; 

X{1} = [randn(15, 1)*0.03 - 1; randn(15, 1)*0.03+1]*10;
X{2} = randn(30, 1)*2;
X{3} = (rand(30, 1)-.5)*30;
noiseLevel = 0.1;
y{1} = cos(pi/2*X{1}/10) + randn(30, 1)*noiseLevel;
y{2} = sin(pi/2*X{2}/10) + randn(30, 1)*noiseLevel;
y{3} = -sin(pi/2*X{3}/10) + randn(30, 1)*noiseLevel;

models = psivm(X, y, kernelType, noiseType, selectionCriterion, d);
models = psivmOptimiseIVM(models, 1);


offset = 0.015*30;
fontName = 'times';
fontSize = 28;

figure(1)
[xvals, yvals] = fplot('cos(pi/2*x/10)', [-15 15]);
a = plot(xvals, yvals);
set(a, 'linewidth', 2)
hold on
figure(2)
[xvals, yvals] = fplot('sin(pi/2*x/10)', [-15 15], 'linewidth', 2);
a = plot(xvals, yvals);
set(a, 'linewidth', 2)
hold on
figure(3)
[xvals, yvals] = fplot('-sin(pi/2*x/10)', [-15 15], 'linewidth', 2);
a = plot(xvals, yvals);
set(a, 'linewidth', 2)
hold on

for taskNo = 1:3

  figure(taskNo)
  cros = plot(X{taskNo}, y{taskNo}, 'rx');
  set(cros, 'linewidth', 3);
  set(cros, 'markersize', 10);
  
  set(gca, 'ylim', [-1.2 1.2], ...
	   'fontSize', fontSize, ...
	   'fontname', fontName, ...
	   'xtick', [-15  -5 0 5 15], ...
	   'ytick', [-1 0 1]);
  
  kVals = find(models.taskOrder == taskNo);
  for i = 1:length(kVals)
    x = X{taskNo}(models.task(taskNo).I(i));
    yval = y{taskNo}(models.task(taskNo).I(i));
    a = plot(x, yval, 'o');
    set(a, 'markersize', 15, 'linewidth', 2)
    set(a, 'erasemode', 'xor')
    %    b = text(x+offset, yval, num2str(kVals(i)));
%    set(b, 'fontName', fontName, 'fontSize', fontSize);
  end
    zeroaxes(gca, 0.0025, 28, 'times')
%    ppos = get(gcf, 'paperposition')
%    ppos(4) = ppos(4)/2;
%    set(gcf, 'paperposition', ppos)
    pos = get(gca, 'position')
    pos(4) = pos(4)/2;
    set(gca, 'position', pos)
    print('-deps', ['../tex/diagrams/sine' num2str(taskNo)])
end