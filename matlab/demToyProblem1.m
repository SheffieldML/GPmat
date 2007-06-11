% DEMTOYPROBLEM1 Generate an artifical data set and solve with GPSIM.

% GPSIM


rand('seed', 1e5);
randn('seed', 1e5);
% protein function is of the form \sum_i \alpha_i
% \exp(-(t-\mu_i)/sigma_i^2)
alpha = [1.5 1 .5]*2;
mu =[0.125 0.8 1.3]*4;
sigma = [0.4 0.3 0.4]*2;

% Properties of genes:
B = [2.6 1.5 0.5 0.2 1.35]; % From Martino paper ... but don't know the scale
S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                             % know the scale, because p21 is
                             % fixed to 1.
D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.
% B = [2 4 .06];
% S = 2*[1 2 1];
% D = [1 1 0.01];
t = linspace(0, 14, 100)';


interSpace = t(2)-t(1);
truef = gpsimArtificialProtein(t, alpha, mu, sigma);
truey = gpsimArtificialGenes(t, alpha, mu, sigma, B, S, D);
datat = [0:1:12]';
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
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity to 1
options.fix(1).index = 8;
options.fix(1).value = expTransform(0.8, 'xtoa');;
options.fix(2).index = 9;
options.fix(2).value = expTransform(1, 'xtoa');;
options.fix(3).index = 15;
options.fix(3).value = expTransform(0.2, 'xtoa');;


model = gpsimCreate(5, 1, datat, dataY, avals./(bvals.*bvals), options);
model = gpsimOptimise(model, 1, 2000);

predt = [linspace(0, 14, 100)]';
proteinKern = kernCreate(model.t, 'rbf');
proteinKern.inverseWidth = ...
    model.kern.comp{1}.inverseWidth;
K = [];
for i=1:model.kern.numBlocks
  K = [K; simXrbfKernCompute(model.kern.comp{i}, proteinKern, ...
                             model.t, predt)];
  
end

obsY = model.y;
startInd = 1;
for i = 1:5
  endInd = i*length(datat);
  obsY(startInd:endInd) = obsY(startInd:endInd)-model.mu(i);
  startInd = endInd + 1;
end

predF = K'*model.invK*obsY;
varF = kernDiagCompute(proteinKern, predt) - sum(K.*(model.invK*K), 1)';

% Take out predictions at data points.
% Use them to get the scale for the other data.
dataF = predF(end-6:end);
dataVarF = varF(end-6:end);
predF(end-6:end) = [];
varF(end-6:end) = [];
predt(end-6:end) = [];
scalePred = sqrt(var(dataF));


figure
lin1 = plot(t, truey);
hold on
lin2 = [plot(datat, dataY, 'x')]
set(lin1, 'lineWidth', 2);
set(lin1, 'markersize', 20);
set(lin2, 'lineWidth', 4);
set(lin2, 'markersize', 20);
set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
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
set(gcf, 'paperposition', origpos)




figure, lin = plot(predt, predF, '-');
hold on,
bh = plot(predt, predF + 2*sqrt(varF), '--');
bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
set(bh, 'lineWidth', 3);
set(lin, 'lineWidth', 4);
set(lin, 'markersize', 20);
set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
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
set(gcf, 'paperposition', origpos)
yLim = get(gca, 'ylim');
  
figure, lin = plot(t, truef, 'r-');
set(lin, 'lineWidth', 4);
set(lin, 'markersize', 20);
set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
set(gca, 'ylim', yLim);
fileName = ['demToyProblem1_true'];
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
set(gcf, 'paperposition', origpos)



save demToyProblem1.mat