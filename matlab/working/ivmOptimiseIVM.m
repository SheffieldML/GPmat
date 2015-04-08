function model = ivmOptimiseIVM(model, display)

% IVMOPTIMISEIVM Optimises an IVM model.

% IVM

if nargin < 2
  display = 1;
end

model = ivmInit(model);
numData = size(model.X, 1);
model.inactiveIndex = 1:numData;

% Set first infoChange to NaN
infoChange(1) = NaN;

for k = 1:model.d
  delta = computeInfoChange(model);
  [infoChange(k), indexSelect] = max(delta);
  if sum(delta==infoChange(k))==length(delta);
    indexSelect = ceil(rand(1)*length(delta));
  end
  i = model.inactiveIndex(indexSelect);
  selectedNum = length(model.activeIndex);
  if selectedNum > 0
    
    switch model.noiseType
      
     case 'probit'
      % update p_i and m_i
      model = probitUpdateSites(model, i);
     case 'gaussian'
      model = gaussianUpdateParams(model, i);
    end
    
    model = updateCholesky(model, i, k, model.inactiveIndex);
    
    switch model.noiseType
     case 'probit'    
      model.h(model.inactiveIndex) = model.h(model.inactiveIndex) + model.alpha(i) * model.L(end,end)*1/ ...
	  sqrt(model.sitePrecision(i))*model.M(k, model.inactiveIndex)';
      model = probitUpdateParams(model, model.inactiveIndex);
     case 'gaussian'
      model = gaussianUpdateParams(model, model.inactiveIndex);
    end
  else
    switch model.noiseType
     case 'probit'
      % update site i's mean and precision
      model = probitUpdateSites(model, i);
     case 'gaussian'
      model = gaussianUpdateParams(model, i);
    end
    
    model = updateCholesky(model, i, k, model.inactiveIndex);
    switch model.noiseType
     case 'probit'
      % Update model
      model.h = model.h + model.alpha(i)*model.L*1/sqrt(model.sitePrecision(i))*model.M';
      model = probitUpdateParams(model);
      
     case 'gaussian'
      model = gaussianUpdateParams(model);
      
    end
  end
  % Remove point from the non-active set and place in the active.
  model.inactiveIndex(indexSelect) = [];
  model.activeIndex = [model.activeIndex; i];
  
  switch model.noiseType
    
   case 'probit'
    logLikelihoods = log(cummGaussian(model.z));
    dLogLikelihood(k) = sum(logLikelihoods);
    logLikelihoodRemain(k) = sum(logLikelihoods(model.inactiveIndex));
    falsePositives(k) = sum(sign(model.h(model.inactiveIndex))~=model.y(model.inactiveIndex) & model.y(model.inactiveIndex)==1);
    trueNegatives(k) = sum(sign(model.h(model.inactiveIndex))~=model.y(model.inactiveIndex) & model.y(model.inactiveIndex)==-1);
   case 'gaussian'
    logLikelihoods = 0;
    dLogLikelihood(k) = 0;
    logLikelihoodRemain(k) = 0;
    falsePositives(k) = 0;
    trueNegatives(k) = 0;
  end
  if display
    switch model.noiseType
     case 'gaussian'
      fprintf('%ith inclusion, remaining log likelihood %2.4f.\n', k, logLikelihoodRemain(k));
     case 'probit'
      fprintf(['%ith inclusion, remaining log Likelihood %2.4f, falsePos %2.4f, trueNeg' ...
	       ' %2.4f\n'], k, logLikelihoodRemain(k), falsePositives(k)./sum(model.y==1), trueNegatives(k)./sum(model.y==-1));
    end
    if display > 1
      if size(model.X, 2) == 2
	figure(2)
	plot(logLikelihoodRemain)
	
	figure(1)
	a = plot(model.X(i, 1), model.X(i, 2), 'o');
	set(a, 'erasemode', 'xor')
	b = text(model.X(i, 1)+0.1, model.X(i, 2), num2str(k));
	set(b, 'erasemode', 'xor')
      else
	subplot(10, 10, rem(k-1, 100)+1);
	image(round(reshape(model.X(i, :), 20, 20)*64))
	axis image
	axis off
      end
      drawnow 
    end
  end
  
  %   if k < initialSelection % Continue selecting from entire data-set
  %   else
  
  %     retainNumber = ceil(model.tau*model.m);
  %     [sortDelta, sortIndex] = sort(delta);
  
%     % Remove low info points
%     model.inactiveIndex = model.inactiveIndex(sortIndex(1:retainNumber));
%     removeIndices = sortIndex(retainNumber+1:end);
    
%     % Put them in R
%     R = [R model.inactiveIndex(sortIndex(retainNumber+1:end))];

%     % Take random sample back from R
%     addNumber = model.m - retainNumber;
%     addIndex = randperm(length(R));
%     addIndex = addIndex(1:addNumber);
%     model.inactiveIndex = [model.inactiveIndex R(addIndex)];

%     % Remove those from R
%     R(addIndex) = [];
    
%     % Now need to update values for the new points which are model.inactiveIndex(retainNumber+1:end)
%     M
%     diagA(removeIndices) = [];
%     h(removeIndices) =[];

  
%   end
end






