% DEMTOYPROBLEM1 Generate an artifical data set and solve with GPSIM.

% GPSIM

maxPredVal = 18;
bw = true;

rand('seed', 1e5);
randn('seed', 1e5);
% protein function is of the form \sum_i \alpha_i
% \exp(-(t-\mu_i)/sigma_i^2)
alpha = [1.5 .8 .8]*2;
mu =[0.125 0.8 1.1]*4+2;
sigma = [0.4 0.5 0.5]*2;

% Properties of genes:
B = [3 3 1 0.2 2.35]; % From Martino paper ... but don't know the scale
S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                             % know the scale, because p21 is
                             % fixed to 1.
D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.
t = linspace(0, maxPredVal, 100)';


interSpace = t(2)-t(1);
truef = gpsimArtificialProtein(t, alpha, mu, sigma);
truey = gpsimArtificialGenes(t, alpha, mu, sigma, B, S, D);
datat = [0:1:16]';
meanDatay = gpsimArtificialGenes(datat, alpha, mu, sigma, ...
                             B, S, D);

bvals = 100;
avals = meanDatay*bvals;

dataY = gammarnd(avals, 1/bvals);


% Check recovery 
fGuess = (diff(truey)/interSpace ...
           +repmat(D, size(truey, 1)-1, 1)...
          .*(truey(1:end-1, :))...
          -repmat(B, size(truey, 1)-1, 1)) ...
          ./repmat(S, size(truey, 1)-1, 1);





% Get the default options structure.
options = gpsimOptions;
options.fix(1).index = 2;
options.fix(1).value = expTransform(1, 'xtoa');;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
 options.fix(2).index = 9;
 options.fix(2).value = expTransform(0.8, 'xtoa');;
% % the corresponding sensitivity to 1
 options.fix(3).index = 10;
 options.fix(3).value = expTransform(1, 'xtoa');;

% Prior information for the TF
options.proteinPrior = [0]';           % Assuming that f=0 at t=0
options.proteinPriorTimes = [0]';
options.optimiser = 'scg';



model = gpsimCreate(5, 1, datat, dataY, avals./(bvals.*bvals), options);
model = gpsimOptimise(model, 1, 2000);

proteinKern = model.kern.comp{1};
for i = 1:model.kern.numBlocks
  predt{i}  = t;
end


K = kernCompute(model.kern, model.timesCell, predt);

obsY = model.m;
predictions =K'*model.invK*obsY;
predF = predictions(1:100, :);
predY = reshape(predictions(101:end), 100, 5);
predY = predY + repmat(model.mu, 100, 1);
variances = kernDiagCompute(model.kern, predt) - sum(K.*(model.invK*K), 1)';
varF = variances(1:100, :);
varY = reshape(variances(101:end), 100, 5);


figure(1)
clf
hold on
for i = 1:5
  stdVals = sqrt(varY(:, i));
  fillColor = [0.7 0.7 0.7];
  fill([predt{1}; predt{1}(end:-1:1)], ...
       [predY(:, i); predY(end:-1:1, i)] ...
       + 2*[stdVals; -stdVals(end:-1:1)], ...
       fillColor,'EdgeColor',fillColor)
end
lin1 = plot(t, truey);
hold on
lin3 =plot(t, predY);
lin2 = [plot(datat, dataY, '.')]
set(lin1, 'lineWidth', 2);
set(lin2, 'markersize', 20);
set(lin3, 'lineWidth', 2);
set(gca, 'fontname', 'times', 'fontsize', 24, 'xlim', xlim, 'ylim', [-2 ...
                    18])
set(gca, 'ytick', [-2:2:18])
zeroAxes(gca);
if bw
  fileName = ['demToyProblem1bw_genes'];
  set(lin1, 'color', [0 0 0])
  set(lin2, 'color', [0 0 0])
  set(lin3, 'color', [0 0 0])
  set(lin1, 'linestyle', '--')
  print('-deps', ['../tex/diagrams/' fileName]);
else
  fileName = ['demToyProblem1_genes'];
  print('-depsc', ['../tex/diagrams/' fileName]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos);
end


figure(2)
clf
stdVals = sqrt(varF);
fillColor = [0.7 0.7 0.7];
fill([predt{1}; predt{1}(end:-1:1)], ...
     [predF; predF(end:-1:1)] ...
     + 2*[stdVals; -stdVals(end:-1:1)], ...
     fillColor,'EdgeColor',fillColor)
hold on
lin = plot(predt{1}, predF, 'b-');
lin2 = plot(t, truef, 'r-');
set(lin, 'lineWidth', 3);
set(lin2, 'lineWidth', 3);
set(lin, 'markersize', 20);
set(gca, 'fontname', 'times', 'fontsize', 24, 'xlim', xlim, 'ylim', [-1 ...
                    7], 'ytick', [-1:7])
zeroAxes(gca);
if bw
  fileName = ['demToyProblem1bw_infered'];
  set(lin, 'color', [0 0 0])
  set(lin2, 'color', [0 0 0], 'lineStyle', '--')
  print('-deps', ['../tex/diagrams/' fileName]);
else
  fileName = ['demToyProblem1_infered'];
  print('-depsc', ['../tex/diagrams/' fileName]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos);
end

bwColormap = [0 0 0; 1 1 1];
colorColormap = [0 0 1; 1 0 0];
figure(5)
h = bar([model.S; S]')
set(gca, 'ylim', [0 3], 'ytick', [0:3])
ylabel('S')

figure(4)
h = [h; bar([model.D; D]')]
set(gca, 'ylim', [0 1], 'ytick', [0:0.5:1])
ylabel('D')

figure(3)
h = [h bar([model.B; B]')]
set(gca, 'ylim', [0 4], 'ytick', [0:4])
ylabel('B')


for i = 3:5
  figure(i)
  set(gca, 'xticklabel', {'gene 1', 'gene 2', 'gene 3', 'gene 4', ['gene ' ...
                      '5']})  
  if bw
    colormap(bwColormap);
    fileName = ['demToyProblem1bw_bar' num2str(i-2)];
    print('-deps', ['../tex/diagrams/' fileName]);
  else
    colormap(colorColormap);
    fileName = ['demToyProblem1_bar' num2str(i-2)];
    print('-depsc', ['../tex/diagrams/' fileName]);
    pos = get(gcf, 'paperposition')
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    print('-dpng', ['../html/' fileName])
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos);
  end

  set(gca, 'fontname', 'times', 'fontsize', 18)
end

save demToyProblem1.mat