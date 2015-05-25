% Set up files for testing noise C++ code


numProcess = 7;
numData = 40;
types = {'gaussian', 'probit', 'mgaussian', 'ngauss', ...
         'ordered', 'ncnm', {'probit', 'gaussian', 'ncnm'}};

for ind=1:length(types)
  if iscell(types{ind})
    % compound noise type
    type = 'cmpnd';
  else
    type=types{ind};
  end
  % for specific cases
  
  noise2 = noiseCreate(types{ind});
  noise2.C = 10; % ordered
  noise2.numProcess=numProcess;
  noise2.numData = numData;
  noise2 = noiseParamInit(noise2);
  % Set the parameters randomly.
  params = noiseExtractParam(noise2);
  params = randn(size(params))./sqrt(randn(size(params)).^2);
  noise2 = noiseExpandParam(noise2, params);
  
  mu = randn(numData, numProcess).*sqrt(1./(randn(numData,numProcess).^2));
  
  varsigma1 = 1./(randn(numData, numProcess).^2);
  varsigma2 = (randn(numData, numProcess).^2);
  varsigmaSwitch = rand(numData, numProcess)> 0.5;
  varsigma = varsigmaSwitch.*varsigma1 + (1-varsigmaSwitch).*varsigma2;
  y = noiseOut(noise2, mu, varsigma);
  
  mu = randn(numData, numProcess).*sqrt(1./(randn(numData,numProcess).^2));
  
  varsigma1 = 1./(randn(numData, numProcess).^2);
  varsigma2 = (randn(numData, numProcess).^2);
  varsigmaSwitch = rand(numData, numProcess)> 0.5;
  varsigma = varsigmaSwitch.*varsigma1 + (1-varsigmaSwitch).*varsigma2;
  
  % Test for missing variables
  if noise2.missing
    index = randperm(size(y, 1)*size(y, 2));
    index = index(1:3);
    y(index) = NaN;
  end
  params = noiseExtractParam(noise2);
  L2 = noiseLogLikelihood(noise2, mu, varsigma, y);
  g2 = noiseGradientParam(noise2, mu, varsigma, y);
  [gmu2, gvs2] = noiseGradVals(noise2, mu, varsigma, y);
    
  save([type 'NoiseTest.mat'], 'noise2', 'params', 'varsigma', 'mu', ...
       'g2', 'gmu2', 'gvs2', 'L2', 'y');
end