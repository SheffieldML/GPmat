function models = psivmOptimiseIVM(models, display)

% PSIVMOPTIMISEIVM Does point selection for a point-set IVM model.

% PSIVM

if nargin < 2
  display = 1;
end

numTasks = length(models.task);
for taskNo = 1:numTasks
  models.task(taskNo) = ivmInit(models.task(taskNo), ceil(models.d*1.2/numTasks));
  numDataPerTask(taskNo) = size(models.task(taskNo).X, 1);
  models.task(taskNo).J = 1:numDataPerTask(taskNo);
end

% Set first infoChange to NaN
trackInfoChange(1) = NaN;

for k = 1:models.d
  
  switch models.selectionCriteria
   case 'random'
    while 1
      taskSelect = ceil(rand(1)*numTasks);
      if length(models.task(taskSelect).J)
	break
      end
    end
    jPos = ceil(rand(1)*length(models.task(taskSelect).J));
    i = models.task(taskSelect).J(jPos);
    trackInfoChange(k) = -.5*log2(1-...
				  models.task(taskSelect).varSigma(models.task(taskSelect).J(jPos))...
				  .*models.task(taskSelect).nu(models.task(taskSelect).J(jPos)));
   case 'entropy'
    for taskNo = 1:numTasks
      delta = -.5*log2(1-models.task(taskNo).varSigma(models.task(taskNo).J).*models.task(taskNo).nu(models.task(taskNo).J));
      if ~isempty(delta)
	[infoChange(taskNo), indexSelect(taskNo)] = max(delta);
	if sum(delta==infoChange(taskNo))==length(delta);
	  indexSelect(taskNo) = ceil(rand(1)*length(delta));
	end
      else
	infoChange(taskNo) = NaN;
	indexSelect(taskNo) = [];
      end
    end
    [trackInfoChange(k) , taskSelect] = max(infoChange);
    jPos = indexSelect(taskSelect);
    i = models.task(taskSelect).J(jPos);
  end
  models.taskOrder(k) = taskSelect;
  selectedNum = length(models.task(taskSelect).I);

  models.task(taskSelect) = ivmAddPoint(models.task(taskSelect), i);

  switch models.noiseType
    case 'probit'
     logLikelihoods = log(cummGaussian(models.task(taskSelect).u));
     dLogLikelihood(k) = sum(logLikelihoods);
     logLikelihoodRemain(k) = sum(logLikelihoods(models.task(taskSelect).J));
     falsePositives(k) = sum(sign(models.task(taskSelect).mu(models.task(taskSelect).J))~=models.task(taskSelect).y(models.task(taskSelect).J) & models.task(taskSelect).y(models.task(taskSelect).J)==1);
     trueNegatives(k) = sum(sign(models.task(taskSelect).mu(models.task(taskSelect).J))~=models.task(taskSelect).y(models.task(taskSelect).J) & models.task(taskSelect).y(models.task(taskSelect).J)==-1);
   case 'gaussian'
    logLikelihoods = -.5*log(2*pi) -.5*((models.task(taskSelect).y-models.task(taskSelect).mu).^2)./(1./models.task(taskSelect).beta ...
						  + models.task(taskSelect).varSigma) - .5*log(1./models.task(taskSelect).beta ...
						  +models.task(taskSelect).varSigma);
    dLogLikelihood(k) = sum(logLikelihoods);
    logLikelihoodRemain(k) = sum(logLikelihoods(models.task(taskSelect).J));
  end

  if display
    if ~rem(k, 10)
      switch models.noiseType
       case 'gaussian'
	fprintf('%ith inclusion, remaining log likelihood %2.4f.\n', k, logLikelihoodRemain(k));
       case 'probit'
	fprintf(['%ith inclusion, remaining log Likelihood %2.4f, falsePos %2.4f, trueNeg' ...
		 ' %2.4f\n'], k, logLikelihoodRemain(k), falsePositives(k)/sum(models.task(taskSelect).y==1), trueNegatives(k)/sum(models.task(taskSelect).y==-1));
      end
    end

    if display > 1
      if size(models.task(taskSelect).X, 2) == 2
	figure(2)
	plot(logLikelihoodRemain)
	
	figure(1)
	plot(models.task(taskSelect).X(i, 1), models.task(taskSelect).X(i, 2), 'o')
	text(models.task(taskSelect).X(i, 1)+0.1, models.task(taskSelect).X(i, 2), num2str(k))
      else
	subplot(10, 10, rem(k-1, 100)+1);
	image(round(reshape(models.task(taskSelect).X(i, :), 20, 20)*64))
	axis image
	axis off
      end
      drawnow 
    end
  end
  
end
models.trackInfoChange = trackInfoChange;




