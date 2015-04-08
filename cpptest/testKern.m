% Set up files for testing kernel C++ code

X = randn(100, 4);
X2 = randn(200, 4);
types = {'bias', 'white', 'rbf', 'lin', {'rbf', 'lin', 'bias', 'white'}};

for ind=1:length(types)
  if iscell(types{ind})
    type = 'cmpnd';
  else
    type = types{ind};
  end
  covGrad = randn(size(X, 1));
  covGrad = covGrad + covGrad';
  kern2 = kernCreate(X, types{ind});
  params = randn(1, kern2.nParams);
  kern2 = kernExpandParam(kern2, params);
  K2 = kernCompute(kern2, X);
  K4 = kernCompute(kern2, X, X2);
  k2 = kernDiagCompute(kern2, X);
  g2 = kernGradient(kern2, X, covGrad);
  G2 = kernGradX(kern2, X, X2);
  GD2 = kernDiagGradX(kern2, X);
  
  save([type 'Test.mat'], 'kern2', 'params', 'K2', 'K4', 'k2', 'X', 'X2', ...
       'g2', 'G2', 'GD2', 'covGrad');
end