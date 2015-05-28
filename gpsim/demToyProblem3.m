% DEMTOYPROBLEM3 Generate an artifical data set and solve with GPSIM.

% SHEFFIELDML

set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextFontSize',20);
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultAxesLineWidth',1.5);

bw = true;
colordef white
rand('seed', 1e6);
randn('seed', 1e6);
% protein function is of the form \sum_i \alpha_i
% \exp(-(t-\mu_i)/sigma_i^2)
alpha = [1.5 1.5 .5 .5];
mu = [4 6 8.5 10.5];
sigma = [1.5 1.5 1.5 1.5];

noiseLevel = 0.01;

% Properties of genes:
B = [eps 0.075 0.0025]; 
S = [1.0 0.400 0.4000]; 
D = [1.0 0.050 0.0010];
t = linspace(0, 18, 100)';
numGenes=length(D);
numObs = 2;
interSpace = t(2)-t(1);
truef = gpsimArtificialProtein(t, alpha, mu, sigma);
truey = gpsimArtificialGenes(t, alpha, mu, sigma, B, S, D);
datat = rand(numObs, 1)*16;%[0:1:12]';
dataY = gpsimArtificialGenes(datat, alpha, mu, sigma, ...
                             B, S, D);
% Get the default options structure.
options = gpsimOptions;
options.includeNoise=1;
% Prior information for the TF
options.proteinPrior = [0]';           % Assuming that f=0 at t=0
options.proteinPriorTimes = [0]';
options.optimiser = 'scg';
options.singleNoise = true;

% Fix kernel values to true values.
options.fix(1).index = 3;
options.fix(1).value = expTransform(D(1), 'xtoa');
options.fix(2).index = 4;
options.fix(2).value = expTransform(S(1)*S(1), 'xtoa');
options.fix(3).index = 5;
options.fix(3).value = expTransform(D(2), 'xtoa');
options.fix(4).index = 6;
options.fix(4).value = expTransform(S(2)*S(2), 'xtoa');
options.fix(5).index = 7;
options.fix(5).value = expTransform(D(3), 'xtoa');
options.fix(6).index = 8;
options.fix(6).value = expTransform(S(3)*S(3), 'xtoa');
options.fix(7).index = 2;              % RBF variance
options.fix(7).value = expTransform(1, 'xtoa');
options.fix(8).index = 9;            
options.fix(8).value = expTransform(noiseLevel*noiseLevel, 'xtoa');
options.fix(9).index = 10;           
options.fix(9).value = expTransform(B(1), 'xtoa');
options.fix(10).index = 11;           
options.fix(10).value = expTransform(B(2), 'xtoa');
options.fix(11).index = 12;           
options.fix(11).value = expTransform(B(3), 'xtoa');


model = gpsimCreate(numGenes, 1, datat, dataY, zeros(size(dataY)), options);
startInd = 1;
endInd = 0;
model.timesCell{2} = [8 16]';
model.timesCell{3} = [6 12]';
model.timesCell{4} = [12 14]';
for i = 2:numGenes+1
  endInd = endInd + numObs;
  %model.timesCell{i} = rand(numObs, 1)*16;
  temp = gpsimArtificialGenes(model.timesCell{i}, alpha, mu, sigma, ...
                             B, S, D);
  model.y(startInd:endInd)=temp(:, i-1)+ randn(numObs, 1)*noiseLevel;
  startInd = endInd + 1;
end

model = gpsimOptimise(model, 1, 100);




K = kernCompute(model.kern, model.timesCell);
invK = pdinv(K);
predt = cell(3, 1);
for i=1:length(model.timesCell)
  predt{i} = t;
end
Kstar = kernCompute(model.kern.comp{1}, predt, model.timesCell);

obsY = model.m;
startInd = 1;
endInd = 0;
for i = 1:length(model.timesCell)
  endInd = endInd+length(model.timesCell{i});
  startInd = endInd + 1;
end

invK = pdinv(K);
ypred = reshape(real(Kstar*invK*obsY), 100, 4);
yvar = reshape(real(kernDiagCompute(model.kern.comp{1}, predt) - sum(Kstar'.* ...
                                                  (invK*Kstar'), 1)'), ...
               100, 4);
KstarStar = kernCompute(model.kern.comp{1}, predt);
ycovar = KstarStar - Kstar*invK*Kstar';
numSamps = 10;
ycovar = 0.5*(ycovar + ycovar');
ysamp = real(gsamp(ypred, ycovar+eye(size(ycovar))*1e-6, numSamps))';
  figure(1)
  clf;
  hold on
startInd = 1;
endInd = 0;
dataPointsHandle = [];
for j = 1:numGenes
  endInd = endInd + numObs;
  dataPointsHandle = [dataPointsHandle plot(model.timesCell{j+1}, model.y(startInd:endInd), '.')];
  startInd = endInd + 1;
end
samplesHandle = [];
trueValueHandle = [];
for j = 1:numGenes
  for i = 1:numSamps
    yvals = reshape(ysamp(:, i), 100, 4);
    samplesHandle = [samplesHandle plot(predt{1}, yvals(:, j+1)+model.mu(j))];
  end
  trueValueHandle = [trueValueHandle plot(predt{1}, truey(:, j))];
end
figure(2)
clf
hold on
for i = 1:numSamps
  yvals = reshape(ysamp(:, i), 100, 4);
  samplesHandle = [samplesHandle plot(predt{1}, yvals(:, 1))];
  trueValueHandle = [trueValueHandle plot(predt{1}, truef)];  
end

for i = 2:length(model.timesCell)
  ypred(:, i) = ypred(:, i) +model.mu(i-1);
end
figure(3)
clf
hold on
fillColor = [0.7 0.7 0.7];
for i = 2:size(ypred, 2)
  fill([predt{1}; predt{1}(end:-1:1)], ...
       [ypred(:, i); ypred(end:-1:1, i)] ...
       + 2*[sqrt(yvar(:, i)); -sqrt(yvar(end:-1:1, i))], ...
       fillColor,'EdgeColor',fillColor)
end
meanHandle = plot(predt{1}, ypred(:, 2:end));
% errorBarHandle = plot(predt{1}, ypred(:, 2:end)+2*sqrt(yvar(:, 2:end)), '--');
% errorBarHandle = [errorBarHandle; plot(predt{1}, ypred(:, 2:end)-2*sqrt(yvar(:, 2:end)), '--')];

startInd = 1;
endInd = 0;
for i = 1:numGenes
  endInd = endInd + numObs;
  dataPointsHandle = [dataPointsHandle plot(model.timesCell{i+1}, model.y(startInd:endInd), '.')];
  startInd=endInd + 1;
end

figure(4)
clf
hold on

fillColor = [0.7 0.7 0.7];
fill([predt{1}; predt{1}(end:-1:1)], ...
     [ypred(:, 1); ypred(end:-1:1, 1)] ...
     + 2*[sqrt(yvar(:, 1)); -sqrt(yvar(end:-1:1, 1))], ...
     fillColor,'EdgeColor',fillColor)
hold on,
meanHandle = [meanHandle; plot(predt{1}, ypred(:, 1))];
%errorBarHandle = [errorBarHandle; plot(predt{1}, ypred(:, 1)+2*sqrt(yvar(:, 1)), '--')];
%errorBarHandle = [errorBarHandle; plot(predt{1}, ypred(:, 1)-2*sqrt(yvar(:, 1)), '--')];

set(dataPointsHandle, 'markersize', 30);
set(trueValueHandle, 'linewidth', 3);
set(samplesHandle, 'linewidth', 1);
%set(errorBarHandle, 'linewidth', 2);
set(meanHandle, 'linewidth', 3);
if bw
  set(dataPointsHandle, 'color', [0 0 0])
  set(trueValueHandle, 'color', [0 0 0])
  set(samplesHandle, 'color', [0 0 0], 'linestyle', '--');
  %set(errorBarHandle, 'color', [0 0 0])
  set(meanHandle, 'color', [0 0 0])
  set(samplesHandle, 'color', [0 0 0], 'linestyle', '--');
else
  set(dataPointsHandle, 'color', [0 0 1])
  set(trueValueHandle, 'color', [0 0 1])
  set(samplesHandle, 'color', [1 0 0]);
end

for i = 1:4
  figure(i);
  set(gca, 'ylim', [-1 7]);
  %zeroAxes(gca);
  if bw
    print('-deps', ['../tex/diagrams/demToyProblem3bw_' num2str(i) ...
                    '.eps']);
  else
    fileName = ['demToyProblem3_' num2str(i)];
    print('-depsc', ['../tex/diagrams/' fileName '.eps']);
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    print('-dpng', ['../html/' fileName])
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos)
  end
  
end



figure(5)
clf
timage = linspace(0, 10, 11)';
counter = 1;
for i = 1:4
  for j = 1:4
    subplot('position', [0.1 0.1 0 0]+ 0.8*[(j-1)*0.25 0.75-(i-1)*0.25 0.23 0.23]) 
    Kimage = multiKernComputeBlock(model.kern.comp{1}, timage, ...
                                   i, j);
    if j>i 
      Kimage = Kimage';
    end
    imagesc(timage, timage, real(Kimage));
    if i~=4
      set(gca, 'xtick', [])
    end
    if j~=1
      set(gca, 'ytick', [])
    end
    %axis equal
    set(gcf, 'paperposition', [0.25 0.25 7 7])
    if bw
      colormap gray
      fileName = ['demToyProblem3bwImage'];
      print('-deps', ['../tex/diagrams/' fileName '.eps']);
    else
      fileName = ['demToyProblem3Image'];
      print('-depsc', ['../tex/diagrams/' fileName '.eps']);
      pos = get(gcf, 'paperposition');
      origpos = pos;
      pos(3) = pos(3)/2;
      pos(4) = pos(4)/2;
      set(gcf, 'paperposition', pos);
      lineWidth = get(gca, 'lineWidth');
      set(gca, 'lineWidth', lineWidth*2);
      print('-dpng', ['../html/' fileName])
      set(gca, 'lineWidth', lineWidth);
      set(gcf, 'paperposition', origpos)
    end
    %colorbar
  end
end

save demToyProblem3.mat
