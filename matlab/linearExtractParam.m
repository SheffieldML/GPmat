function [params, names] = linearExtractParam(model,dim);

% LINEAREXTRACTPARAM Extract weights from a linear model.
%
% MODIFICATIONS : Carl Henrik Ek, 2009
%
% MLTOOLS

if(nargin<2)
  params = [model.W(:)' model.b];
else
  params = model.W(:,dim);
  params = [params(:)' model.b(dim)];
end

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
  
