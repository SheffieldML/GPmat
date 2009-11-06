% SIMSAMPLE Sample from SIM kernel

% KERN

colordef white
randn('seed', 1e8)
rand('seed', 1e8)

numSamp = 100;
t = linspace(0, 5, numSamp)';
kern = kernCreate(t, {'multi', 'rbf', 'sim', 'sim', 'sim'});
for i = 1:length(kern.comp)
  kern.comp{i}.inverseWidth = 1/(0.75*0.75);
end
kern.comp{2}.decay = 5;
kern.comp{2}.variance = 25;
kern.comp{3}.decay = 1;
kern.comp{3}.variance = 1;
kern.comp{4}.decay = 0.5;
kern.comp{4}.variance = 0.25;

params = kernExtractParam(kern);
kern = kernExpandParam(kern, params);

K = kernCompute(kern, t);

imagesc(K, [-1.1 1.1]);
handle = [];
fontSize = 24;
set(gca, 'fontname', 'times')
set(gca, 'fontsize', fontSize)
handle = [handle; text('Interpreter', 'latex', 'string', '$$f(t)$$', 'position', ...
            [50 410], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_1(t)$$', 'position', ...
            [150 410], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_2(t)$$', 'position', ...
            [250 410], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_3(t)$$', 'position', ...
            [350 410], 'fontsize', fontSize, 'horizontalalignment', 'center')];
set(handle, 'horizontalalignment', 'right')
handle = [handle; text('Interpreter', 'latex', 'string', '$$f(t)$$', 'position', ...
            [-10 50], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_1(t)$$', 'position', ...
            [-10 150], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_2(t)$$', 'position', ...
            [-10 250], 'fontsize', fontSize, 'horizontalalignment', 'center')];
handle = [handle; text('Interpreter', 'latex', 'string', '$$x_3(t)$$', 'position', ...
            [-10 350], 'fontsize', fontSize, 'horizontalalignment', 'center')];
axis off
colorbar
if exist('printDiagram') && printDiagram
  printPlot('simTestKernelImage', '../tex/diagrams', '../html');
end

numSamp2 = 150;
y = gsamp(zeros(size(K), 1), 0.5*(real(K) + real(K')), numSamp2)';
t = t - 0;
counter = 0;
for i = 1:numSamp2
  x = reshape(y(:, i), numSamp, length(kern.comp));
  if any(x(:, 1)<0)
    continue
  end
  figure
  counter = counter + 1;
  a = plot(t, x(:, 1), 'k-');
  hold on
  a = [a plot(t, x(:, 2), 'r-')];
  a = [a plot(t, x(:, 3), 'g-')];
  a = [a plot(t, x(:, 4), 'b-')];
  set(a, 'linewidth', 2);
  set(gca, 'fontname', 'arial');
  set(gca, 'fontsize', 14);
  if exist('printDiagram') && printDiagram
    printPlot(['simSample'  num2str(counter)], '../tex/diagrams', '../html');
  end
end
  
