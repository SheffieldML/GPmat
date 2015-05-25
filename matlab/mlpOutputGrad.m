function g = mlpOutputGrad(model, X)

% MLPOUTPUTGRAD Evaluate derivatives of mlp model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of a multi-layer perceptron's
% outputs with respect to the parameters of the multi-layer
% perceptron. Currently it simply wraps the NETLAB mlpderiv
% function.
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the multi-layer
% perceptron with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : mlpCreate, mlpderiv
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

if length(model.hiddenDim) == 1
  g = mlpderiv(model, X);
else
  
  [Y, G, A] = mlpOut(model, X);
  gw = cell(1, length(model.w));
  for i = 1:length(gw)
    gw{i} = zeros(size(model.w{i}));
  end   
  for i = 1:length(G);
    WdG{i} = (1-G{i}.*G{i})*w{i};
  end
  for k = 1:model.outputDim
    gw{end}(:, k) = Z{end}(:, k);
    gb{end} = 1;
    for i = length(model.w)-1:-1:1
      %gw{i} = 
      error('Not yet implemented');
    end
  end

% Evaluate second-layer gradients.
gw2 = z'*deltas;
gb2 = sum(deltas, 1);

% Now do the backpropagation.
delhid = deltas*net.w2';
delhid = delhid.*(1.0 - z.*z);

% Finally, evaluate the first-layer gradients.
gw1 = x'*delhid;
gb1 = sum(delhid, 1);

g = [gw1(:)', gb1, gw2(:)', gb2];
end
