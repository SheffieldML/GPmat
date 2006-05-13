function [params, names] = mlpExtractParam(model)

% MLPEXTRACTPARAMS Wrapper for NETLAB's mlppak.

% MLTOOLS

params = mlppak(model);

if nargout > 1
  counter = 0;
  for j = 1:size(model.w1, 2)
    for i = 1:size(model.w1, 1)
      counter = counter + 1;
      names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.b1, 2)
    counter = counter + 1;
    names{counter} = ['Hidden node bias ' num2str(j)];
  end
  for j = 1:size(model.w2, 2)
    for i = 1:size(model.w2, 1)
      counter = counter + 1;
      names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.b2, 2)
    counter = counter + 1;
    names{counter} = ['Output node bias ' num2str(j)];
  end
end