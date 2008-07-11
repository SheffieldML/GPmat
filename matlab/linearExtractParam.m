function [params, names] = linearExtractParam(model);

% LINEAREXTRACTPARAM Extract weights from a linear model.

% MLTOOLS

  params = [model.W(:)' model.b];
if nargout > 1
  counter = 0;
  for j = 1:size(model.W, 2)
    for i = 1:size(model.W, 1)
      counter = counter + 1;
      names{counter} = ['Weight ' num2str(i) '-' num2str(j)];
    end
  end
    for j = 1:size(model.b, 2)
    counter = counter + 1;
    names{counter} = ['Bias ' num2str(j)];
  end
end
  