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
  J{taskNo} = 1:numDataPerTask(taskNo);
end

% Set first infoChange to NaN
trackInfoChange(1) = NaN;

for k = 1:models.d
  
  for taskNo = 1:numTasks
    switch models.selectionCriteria
     case 'entropy'
      delta = -.5*log2(1-models.task(taskNo).diagA(J{taskNo}).*models.task(taskNo).nu(J{taskNo}));
    end
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
  i = J{taskSelect}(indexSelect(taskSelect));
  selectedNum = length(models.task(taskSelect).activeIndex);
  if selectedNum > 0

    switch models.noiseType
     case 'probit'
      % update p_i and m_i
      models.task(taskSelect) = probitUpdateSites(models.task(taskSelect), i);
     case 'gaussian'
      models.task(taskSelect) = gaussianUpdateParams(models.task(taskSelect), i);
    end

    models.task(taskSelect) = updateCholesky(models.task(taskSelect), i,  length(models.task(taskSelect).activeIndex)+1, J{taskSelect});
    
    switch models.noiseType
     case 'probit'    
      models.task(taskSelect).h(J{taskSelect}) = models.task(taskSelect).h(J{taskSelect}) + models.task(taskSelect).alpha(i) * models.task(taskSelect).L(end,end)*1/ ...
	  sqrt(models.task(taskSelect).sitePrecision(i))*models.task(taskSelect).M( length(models.task(taskSelect).activeIndex)+1, J{taskSelect})';
      models.task(taskSelect) = probitUpdateParams(models.task(taskSelect), J{taskSelect});
     case 'gaussian'
      models.task(taskSelect) = gaussianUpdateParams(models.task(taskSelect), J{taskSelect});
    end
  else
    switch models.noiseType
     case 'probit'
      % update site i's mean and precision
      models.task(taskSelect) = probitUpdateSites(models.task(taskSelect), i);
     case 'gaussian'
      models.task(taskSelect) = gaussianUpdateParams(models.task(taskSelect), i);
    end
    
    models.task(taskSelect) = updateCholesky(models.task(taskSelect), i, length(models.task(taskSelect).activeIndex)+1, J{taskSelect});
    
    switch models.noiseType
     case 'probit'
      % Update model

      models.task(taskSelect).h = models.task(taskSelect).h + models.task(taskSelect).alpha(i)*models.task(taskSelect).L*1/sqrt(models.task(taskSelect).sitePrecision(i))*models.task(taskSelect).M';
      models.task(taskSelect) = probitUpdateParams(models.task(taskSelect));
      
     case 'gaussian'
      models.task(taskSelect) = gaussianUpdateParams(models.task(taskSelect));
      
    end
  end
  % Remove point from the non-active set and place in the active.
  J{taskSelect}(indexSelect(taskSelect)) = [];
  models.task(taskSelect).activeIndex = [models.task(taskSelect).activeIndex; i];

  switch models.noiseType
    case 'probit'
     logLikelihoods = log(cummGaussian(models.task(taskSelect).z));
     dLogLikelihood(k) = sum(logLikelihoods);
     logLikelihoodRemain(k) = sum(logLikelihoods(J{taskSelect}));
     falsePositives(k) = sum(sign(models.task(taskSelect).h(J{taskSelect}))~=models.task(taskSelect).y(J{taskSelect}) & models.task(taskSelect).y(J{taskSelect})==1);
     trueNegatives(k) = sum(sign(models.task(taskSelect).h(J{taskSelect}))~=models.task(taskSelect).y(J{taskSelect}) & models.task(taskSelect).y(J{taskSelect})==-1);
   case 'gaussian'
    logLikelihoods = 0;
    dLogLikelihood(k) = 0;
    logLikelihoodRemain(k) = 0;
    falsePositives(k) = 0;
    trueNegatives(k) = 0;
  end
  if display
    if ~rem(k, 10)
      switch models.noiseType
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






