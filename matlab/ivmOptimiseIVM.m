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

% If Randomised greedy selection is being used, use dprime.
if isfield(model, 'dprime')
  dVal = model.dprime;
else
  dVal = model.d;
end

for k = 1:dVal
    
  [indexSelect, infoChange(k)] = selectPoint(model);
  dataIndexSelect = model.J(indexSelect);
  
  model = ivmAddPoint(model, dataIndexSelect);
  
  if display
    logLikelihoods = log(ivmLikelihoods(model));
    dLogLikelihood(k) = sum(sum(logLikelihoods));
    logLikelihoodRemain(k) = sum(sum(logLikelihoods(model.J, :)));
    fprintf('%ith inclusion, remaining log Likelihood %2.4f', ...
            k, logLikelihoodRemain(k))
    switch model.noise.type
     case {'probit', 'heaviside'}
      falsePositives(k) = 0;
      truePositives(k) = 0;
      for i = 1:size(model.y, 2)
        falsePositives(k) = falsePositives(k) ...
            + sum(...
                sign(model.mu(model.J, i)+model.noise.bias(i)) ...
                ~=model.y(model.J, i) & model.y(model.J, i)==-1);
        truePositives(k) = truePositives(k) ...
            + sum(...
                sign(model.mu(model.J, i)+model.noise.bias(i)) ...
                ~=model.y(model.J, i) & model.y(model.J, i)==1);
      end
      fprintf(', falsePos %2.4f, truePos %2.4f\n', ...
              sum(falsePositives(k))./sum(sum(model.y==-1)), ...
              sum(truePositives(k))./sum(sum(model.y==1)));
     otherwise
      fprintf('\n');
    end
    if display > 1
      if size(model.X, 2) == 2
%/~	figure(2)
%	plot(logLikelihoodRemain)
%~/	
	figure(1)
	a = plot(model.X(dataIndexSelect, 1), ...
                 model.X(dataIndexSelect, 2), 'o');
	set(a, 'erasemode', 'xor')
	xlim = get(gca, 'xlim');
	labelGap = (xlim(2) - xlim(1)) * 0.025;
	b = text(model.X(dataIndexSelect, 1)+labelGap, ...
                 model.X(dataIndexSelect, 2), num2str(k));
	set(b, 'erasemode', 'xor')
      else
	subplot(10, 10, rem(k-1, 100)+1);
	image(round(reshape(model.X(dataIndexSelect, :), 20, 20)*64))
	axis image
	axis off
      end
      drawnow 
    end
  end
end

