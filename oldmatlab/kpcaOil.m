% KPCAOIL Demonstrate kernel parameter selection with the sub-set of the oil data.

load ..\3Class

latentDim = 2;

% Select data
numData = 100;
indices = randperm(1000);
indices = indices(1:numData);
Y = DataTrn(indices, :);
lbls = DataTrnLbls(indices, :);


meanData = mean(Y);
Y = Y - repmat(meanData, numData, 1);
beta = 500;
dataDim = size(Y, 2);
for gamma = [0.001 0.01 0.1]
  theta = [gamma 1 beta];
  Ky = computeKernel(Y, theta);
  [U, D] = eig(Ky);
  [Lambda, index] = sort(-diag(D));
  Lambda = -Lambda;
  U = U(:, index);
  
  X = zeros(numData, latentDim);
  L(1, 1) = sqrt(Lambda(1) -1/beta);
  L(2, 2) = sqrt(Lambda(2) -1/beta);
  X(:, 1) = U(:, 1);
  X(:, 2) = U(:, 2);
  X = X*L;
  disp(L)
  invKx = pdinv(X*X'+1/beta*eye(numData));
  
  KyinvKx = Ky*invKx;
  val = -dataDim*numData/2 - dataDim/2*log(det(KyinvKx)) + dataDim/2*trace(KyinvKx);
  returnVal = [];
  
  symbol{1} = 'r+';
  symbol{2} = 'bo';
  symbol{3} = 'mx';
  figure, hold on
  for i = 1:size(X, 1)
    labelNo = find(lbls(i, :));
    plot(X(i, 1), X(i, 2), symbol{labelNo})
  end
  title(num2str(val));
end