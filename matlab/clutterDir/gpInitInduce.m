function X_u = gpInitInduce(model)

% GPINITINDUCE Initialise FGPLVM inducing variables


x = model.X_u(:)';
options = optOptions;
options(9) = 1;
fhandle = str2func(model.optimiser);
x = fhandle('induceObjective', x, options, 'induceGradients', model);
X_u = reshape(x, model.k, model.q);








