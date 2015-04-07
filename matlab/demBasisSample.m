% DEMBASISSAMPLE Do a simple demo of a basis function.

% GPSIM

lengthScale = 1;
locationParams = [2 4 6];
numPoints = 100;
numFuncs = 3;


fontName = 'vera';
randn('seed', 1e6)
rand('seed', 1e6)

minLim = min(locationParams)-2*lengthScale;
maxLim = max(locationParams)+2*lengthScale;

colors = {'r', 'g', 'b'};

figure
hold on
basisFunctions = zeros(numPoints, length(locationParams));
t = linspace(minLim, maxLim, numPoints)';
for i = 1:length(locationParams)
  tcentred = t - locationParams(i);
  basisFunctions(:, i) = exp(-tcentred.*tcentred/((lengthScale^2)))/(sqrt(2*pi)*lengthScale);
  a = plot(t, basisFunctions(:, i), colors{i});
  set(a, 'linewidth', 2);  
end
pos = get(gca, 'position');
origpos = pos;
%pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gca, 'position', pos);
%set(gca, 'fontname', fontName);
set(gca, 'fontsize', 20);
printPlot(['demBasisSample', num2str(0)], '../tex/diagrams/', '../html/')



decays = [0.01 0.1 1]
sensitivity = [1 1 1]
mrnaBasisFunctions = zeros(numPoints, length(decays), length(locationParams));
sigma2 = lengthScale*lengthScale;
for i = 1:length(decays)
  figure
  hold on
  for j = 1:length(locationParams);
    decays2=decays(i)*decays(i);
    tcentred = t - locationParams(j);
    mrnaBasisFunctions(:, i, j) = exp(...
        (decays2*sigma2/4 - decays(i)*(t-locationParams(j))) ...
        +lnDiffCumGaussian(2/sqrt(2)*repmat( ...
            (decays(i)*sigma2 + locationParams(j))/lengthScale, ...
                           size(t, 1), size(t, 2)),...
        -2/sqrt(2)*(t - decays(i)*sigma2 - locationParams(j))/lengthScale));
    a = plot(t, mrnaBasisFunctions(:, i, j), colors{j});
    set(a, 'linewidth', 2);  
  end
  pos = get(gca, 'position');
  origpos = pos;
  %pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gca, 'position', pos);
  %set(gca, 'fontname', fontName);
  set(gca, 'fontsize', 20);
  printPlot(['demBasisSample', num2str(0), '_', num2str(i)], '../tex/diagrams/', '../html/')
end
  
w = randn(length(locationParams), numFuncs);
for i=1:numFuncs
  figure
  f = basisFunctions*w(:, i);
  a = plot(t, f, 'b');
  set(a, 'linewidth', 2);  
  pos = get(gca, 'position');
  origpos = pos;
  pos(4) = pos(4)/2;
  set(gca, 'position', pos);
  disp(w(:, i))
  set(a, 'linewidth', 2);
  %set(gca, 'fontname', fontName);
  set(gca, 'fontsize', 20);
  printPlot(['demBasisSample', num2str(i)], '../tex/diagrams/', '../html/')
end


for i=1:length(decays)
  for j=1:numFuncs
    figure
    f = squeeze(mrnaBasisFunctions(:, i, :))*w(:, j);
    a = plot(t, f, 'b');
    set(a, 'linewidth', 2);  
    pos = get(gca, 'position');
    origpos = pos;
    pos(4) = pos(4)/2;
    set(gca, 'position', pos);
    set(a, 'linewidth', 2);
    %set(gca, 'fontname', fontName);
    set(gca, 'fontsize', 20);
    printPlot(['demBasisSample', num2str(j), '_', num2str(i)], '../tex/diagrams/', '../html/')
  end
end
