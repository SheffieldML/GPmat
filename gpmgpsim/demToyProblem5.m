% DEMTOYPROBLEM5 Generate an artifical data set and solve with GPSIM.

% SHEFFIELDML

load demToyProblem5.mat
bw = true;

% protein function is of the form \sum_i \alpha_i
% \exp(-(t-\mu_i)/sigma_i^2)
alpha = [1.5 1.5 .5 .5];
mu = [4 6 8.5 10.5];
sigma = [1.5 1.5 1.5 1.5];

% Properties of genes:
B = [eps 0.075 0.0025]; 
S = [1.0 0.400 0.4000]; 
D = [1.0 0.050 0.0010];
t = linspace(0, 18, 100)';
numGenes=length(D);
truef = gpsimArtificialProtein(t, alpha, mu, sigma);
truey = gpsimArtificialGenes(t, alpha, mu, sigma, B, S, D);

% Get the default options structure.
options = gpsimOptions;

% Learn noise level
options.includeNoise = true;
options.singleNoise = true;

% Prior information for the TF
options.proteinPrior = [0]';           % Assuming that f=0 at t=0
options.proteinPriorTimes = [0]';
options.optimiser = 'scg';

% Fix RBF variance.
options.fix(1).index = 2;              % RBF variance
options.fix(1).value = expTransform(1, 'xtoa');
for repNum = 10
  counter1 = 0;
  rand('seed', repNum*1e6);
  randn('seed', repNum*1e6);
  for numObs = [2 5 10 15 20]
    counter1 = counter1 + 1;
    counter2 = 0;
    for noiseLevel = [0.0001 0.001 0.01 0.1 1]
      counter2 = counter2+1;
      datat = rand(numObs, 1)*16;
      dataY = gpsimArtificialGenes(datat, alpha, mu, sigma, ...
                                   B, S, D);
      model{counter1, counter2, repNum} = gpsimCreate(numGenes, 1, datat, dataY, zeros(size(dataY)), options);
      startInd = 1;
      endInd = 0;
      for i = 2:numGenes+1
        endInd = endInd + numObs;
        model{counter1, counter2, repNum}.timesCell{i} = rand(numObs, 1)*16;
        temp = gpsimArtificialGenes(model{counter1, counter2, repNum}.timesCell{i}, alpha, mu, sigma, ...
                                    B, S, D);
        model{counter1, counter2, repNum}.y(startInd:endInd)=temp(:, i-1)+ randn(numObs, 1)*sqrt(noiseLevel);
        startInd = endInd + 1;
      end
      
      model{counter1, counter2, repNum} = gpsimOptimise(model{counter1, counter2, repNum}, 1, 5000);
      
      
      
      K = kernCompute(model{counter1, counter2, repNum}.kern, model{counter1, counter2, repNum}.timesCell);
      invK = pdinv(K);
      predt = cell(3, 1);
      for i=1:length(model{counter1, counter2, repNum}.timesCell)
        predt{i} = t;
      end
      Kstar = kernCompute(model{counter1, counter2, repNum}.kern.comp{1}, predt, model{counter1, counter2, repNum}.timesCell);
      
      obsY = model{counter1, counter2, repNum}.m;
      startInd = 1;
      endInd = 0;
      for i = 1:length(model{counter1, counter2, repNum}.timesCell)
        endInd = endInd+length(model{counter1, counter2, repNum}.timesCell{i});
        startInd = endInd + 1;
      end
      
      invK = pdinv(K);
      ypred = reshape(real(Kstar*invK*obsY), 100, 4);
      yvar = reshape(real(kernDiagCompute(model{counter1, counter2, repNum}.kern.comp{1}, predt) - sum(Kstar'.* ...
                                                        (invK*Kstar'), 1)'), ...
                   100, 4);
      KstarStar = kernCompute(model{counter1, counter2, repNum}.kern.comp{1}, predt);
      ycovar = KstarStar - Kstar*invK*Kstar';
      numSamps = 10;
      ycovar = 0.5*(ycovar + ycovar');
      fhat = ypred(:, 1)-truef;
      mse(counter1, counter2, repNum) = sum(fhat.*fhat)/length(fhat);
      [invCov, U] = pdinv(ycovar(1:length(t), 1:length(t)));
      tll(counter1, counter2, repNum) = -length(truef)*.5*log(2*pi)-0.5*logdet(ycovar(1:length(t), 1:length(t)),U) ...
          -0.5*(fhat'*invCov*fhat);
      noiseLevels(counter1, counter2, repNum) = noiseLevel;
      numObservations(counter1, counter2, repNum) = numObs;
      save demToyProblem5.mat model tll mse noiseLevels numObservations
    end
  end
end

mVal= mean(mse(:, :, 1:8), 3)';
stdVal = std(mse(:, :, 1:8),[], 3)';
errorbar(mVal, 2*stdVal)
