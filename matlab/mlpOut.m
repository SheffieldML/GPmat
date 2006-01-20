function Y = mlpOut(model, X);

% MLPOUT Output of an MLP model (wrapper for the NETLAB function mlpfwd).

% MLTOOLS

Y = mlpfwd(model, X);