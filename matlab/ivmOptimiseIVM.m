function model = ivmOptimiseIVM(model, display)

% IVMOPTIMISEIVM Optimises an IVM model.

% IVM

if nargin < 2
  display = 1;
end

model = ivmInit(model);
numData = size(model.X, 1);
model.J = 1:numData;

% Set first infoChange to NaN
infoChange(1) = NaN;

% If RGS is being used, use dprime.
if isfield(model, 'dprime')
  dVal = model.dprime;
else
  dVal = model.d;
end

for k = 1:dVal
    
  [indexSelect, infoChange(k)] = selectPoint(model);
  i = model.J(indexSelect);
  
  model = ivmAddPoint(model, i);
  
  switch model.noise.type
   
   case 'probit'
     logLikelihoods = log(cumGaussian(model.u));
     dLogLikelihood(k) = sum(sum(logLikelihoods));
     logLikelihoodRemain(k) = sum(sum(logLikelihoods(model.J, :)));
     falsePositives(k, :) = sum(sign(model.mu(model.J, :)+model.noise.bias)~=model.y(model.J, :) ...
			     & model.y(model.J, :)==1);
     trueNegatives(k, :) = sum(sign(model.mu(model.J, :)+model.noise.bias)~=model.y(model.J, :) ...
			       & model.y(model.J, :)==-1);
   case 'heaviside'
     logLikelihoods = log(cumGaussian(model.u)*(1-2*model.noise.eta)+model.noise.eta);
     dLogLikelihood(k) = sum(sum(logLikelihoods));
     logLikelihoodRemain(k) = sum(sum(logLikelihoods(model.J, :)));
     falsePositives(k, :) = sum(sign(model.mu(model.J, :)+model.noise.bias)~=model.y(model.J, :) & model.y(model.J, :)==1);
     trueNegatives(k, :) = sum(sign(model.mu(model.J, :)+model.noise.bias)~=model.y(model.J, :) & model.y(model.J, :)==-1);
   case 'ordered'
    % Not yet implemented
    logLikelihoods = zeros(size(model.y));
    logLikelihoodRemain(k) = sum(sum(logLikelihoods(model.J, :)));
   case 'gaussian'
    logLikelihoods = -.5*log(2*pi) -.5*((model.y-model.mu).^2)./(model.noise.sigma2 ...
						  + model.varSigma) - .5*log(model.noise.sigma2 ...
						  +model.varSigma);
    dLogLikelihood(k) = sum(sum(logLikelihoods));
    logLikelihoodRemain(k) = sum(sum(logLikelihoods(model.J, :)));
  end
  if display
    switch model.noise.type
     case 'gaussian'
      fprintf('%ith inclusion, remaining log likelihood %2.4f.\n', k, logLikelihoodRemain(k));
     case {'probit', 'heaviside'}
      fprintf(['%ith inclusion, remaining log Likelihood %2.4f, falsePos %2.4f, trueNeg' ...
	       ' %2.4f\n'], k, logLikelihoodRemain(k), sum(falsePositives(k))./sum(sum(model.y==1)), sum(trueNegatives(k))./sum(sum(model.y==-1)));
    end
    if display > 1
      if size(model.X, 2) == 2
	figure(2)
	plot(logLikelihoodRemain)
	
	figure(1)
	a = plot(model.X(i, 1), model.X(i, 2), 'o');
	set(a, 'erasemode', 'xor')
	xlim = get(gca, 'xlim');
	labelGap = (xlim(2) - xlim(1)) * 0.025;
	b = text(model.X(i, 1)+labelGap, model.X(i, 2), num2str(k));
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
end

