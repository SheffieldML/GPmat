% DEMRATEMAPS2
% Fix seeds
importTool('speech')

randn('seed', 1e5);
rand('seed', 1e5);

display = 0;

dataSetName = 'ratemaps';
experimentNo = 2;


baseDir = datasetsDirectory;
dirSep = filesep;
t = linspace(-pi, pi, 256);
window = (cos(t) + 1)/2;
[x, fs] = wavread([baseDir 'osr' dirSep 'OSR_uk_000_0020_8k.wav']);
[r, theta] = stft(x, fs, 256, window);
retainFreq = 128
Y = [r(:, 1:retainFreq) theta(:, 1:retainFreq)];
seq = [80, 290 
       340, 520 
       595, 750
       810, 980
       1025, 1200
       1275, 1430
       1510 1680
       1760 1950
       2030 2200
       2240 2450];
figure(1)
for i = 1:size(seq, 1)
  subplot(size(seq, 1), 1, i)
  pcolor(seq(i, 1):seq(i, 2), 1:retainFreq, log10(Y(seq(i, 1):seq(i, 2), ...
                                                    1:retainFreq )'))
  shading interp
end
for i = 1; %:size(seq, 1)
  rReconstruct = Y(seq(i, 1):seq(i, 2), [1:retainFreq retainFreq:-1:1]);
  thetaReconstruct = [Y(seq(i, 1):seq(i, 2), retainFreq+1:2*retainFreq) ...
                      -Y(seq(i, 1):seq(i, 2), [retainFreq+1 retainFreq*2:-1:retainFreq+2])];
  
  xReconstruct = istft(rReconstruct, ...
                       thetaReconstruct,...
                       fs, ...
                     length(window), window);
  wavplay(xReconstruct, 8000)
end
%seq
t = [];
for i = 1:3
  t = [t; [seq(i, 1):seq(i, 2)]'];
end

Y = log10(Y(t, 1:retainFreq));


latentDim = 4;
options = fgplvmOptions('fitc');
model = fgplvmCreate(latentDim, size(Y, 2), Y, options);
model.scale = sqrt(var(Y));
model.m = gpComputeM(model);
optionsDyn = gpOptions('fitc');
optionsDyn.kern = kernCreate(model.X, {'rbf', 'white'});
optionsDyn.kern.comp{1}.inverseWidth = 0.04;
% This gives signal to noise of 0.1:1e-3 or 100:1.
optionsDyn.kern.comp{1}.variance = 1;
optionsDyn.kern.comp{2}.variance = 1e-3^2;
diff = 1;
learn = 0;
model = fgplvmAddDynamics(model, 'gpTime', optionsDyn, t);
model.dynamics.beta = 1e4;
model.optimiser = 'conjgrad';

% Fit the GP latent variable model
model = fgplvmOptimise(model, 1, 1000);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
%save(['dem' capName num2str(experimentNo) '.mat'], 'X', 'kern', 'noise', 'ivmInfo');
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

%demRateMaps1Project
% Load the results and display dynamically.
% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'spectrum')
%fgplvmResultsDynamic(dataSetName, experimentNo, 'spectrum', 'prepSpectrum', ...
 %                   500);


% Load the results and display statically.
% gplvmResultsStatic(dataSetName, experimentNo, 'vector');

% Load the results and display as scatter plot
% gplvmResultsStatic(dataSetName, experimentNo, 'none')

