function ivmSelectVisualise(model, display, k, dataIndexSelect)

% IVMSELECTVISUALISE Visualise the selected point.
% FORMAT
% DESC visualises the IVM learning by showing which point has been
% selected.
% ARG model : the model for which the visualisation is occuring.
% ARG display : the verbosity level (whether or not to display in a
% plot, 1 - don't plot, 2 - attempt to plot).
% ARG k : the point inclusion number that we are displaying.
% ARG index : the index of the data point that is being included.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : ivmCreate, ivmLogLikelihood, ivmOptions

% IVM

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
    type = model.noise.type;
    if strcmp(type, 'cmpnd')
      type = model.noise.comp{1}.type;
    end
    switch type
     case {'probit', 'ordered', 'ncnm'}
      figure(1)
      a = plot(model.X(dataIndexSelect, 1), ...
               model.X(dataIndexSelect, 2), '.', ...
               'markersize', 20);
      if ~isoctave
        set(a, 'eraseMode', 'xor');
      end
      drawnow
     case {'gaussian', 'ngauss'}
      figure(1)
        a = plot3(model.X(dataIndexSelect, 1), ...
                model.X(dataIndexSelect, 2), model.y(dataIndexSelect), ...
                '.', 'markersize', 20);
        if ~isoctave
          set(a, 'erasemode', 'xor');
        end
      %      xlim = get(gca, 'xlim');
      %      labelGap = (xlim(2) - xlim(1)) * 0.025;
      %      b = text(model.X(dataIndexSelect, 1)+labelGap, ...
      %               model.X(dataIndexSelect, 2), model.y(dataIndexSelect), ...
      %               num2str(k), 'erasemode', 'xor');
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
