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
    
  [indexSelect, infoChange(k)] = ivmSelectPoint(model);
  dataIndexSelect = model.J(indexSelect);
  
  model = ivmAddPoint(model, dataIndexSelect);
  if display
    logLikelihood = ivmLogLikelihood(model);
    fprintf('%ith inclusion, remaining log Likelihood %2.4f', ...
            k, logLikelihood)
    switch model.noise.type
     case {'probit'}
      falsePositives(k) = 0;
      truePositives(k) = 0;
      for i = 1:size(model.y, 2)
        falsePositives(k) = falsePositives(k) ...
            + sum(...
                sign(model.mu(:, i)+model.noise.bias(i)) ...
                ~=model.y(:, i) & model.y(:, i)==-1);
        truePositives(k) = truePositives(k) ...
            + sum(...
                sign(model.mu(:, i)+model.noise.bias(i)) ...
                ~=model.y(:, i) & model.y(:, i)==1);
      end
      fprintf(', falsePos %2.4f, truePos %2.4f\n', ...
              sum(falsePositives(k))./sum(sum(model.y==-1)), ...
              sum(truePositives(k))./sum(sum(model.y==1)));
     otherwise
      fprintf('\n');
    end
    if display > 1
      if size(model.X, 2) == 2
        switch model.noise.type
          case {'probit', 'ordered'}
           figure(1)
           a = plot(model.X(dataIndexSelect, 1), ...
                    model.X(dataIndexSelect, 2), 'o', 'erasemode', 'xor');
           xlim = get(gca, 'xlim');
           labelGap = (xlim(2) - xlim(1)) * 0.025;
           b = text(model.X(dataIndexSelect, 1)+labelGap, ...
                    model.X(dataIndexSelect, 2), num2str(k), 'erasemode', ...
                    'xor');
           drawnow
         case 'gaussian'
           figure(1)
           a = plot3(model.X(dataIndexSelect, 1), ...
                    model.X(dataIndexSelect, 2), model.y(dataIndexSelect), ...
                     'o', 'erasemode', 'xor');
           xlim = get(gca, 'xlim');
           labelGap = (xlim(2) - xlim(1)) * 0.025;
           b = text(model.X(dataIndexSelect, 1)+labelGap, ...
                    model.X(dataIndexSelect, 2), model.y(dataIndexSelect), ...
                    num2str(k), 'erasemode', 'xor');
           drawnow
        end
      else
	subplot(10, 10, rem(k-1, 100)+1);
	image(round(reshape(model.X(dataIndexSelect, :), 20, 20)*64))
	axis image
	axis off
        drawnow
      end
      
    end
  end
end

