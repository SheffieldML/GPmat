function model = orderedUpdateParams(model, index)

% ORDEREDUPDATEPARAMS Update parameters for ordered categorical noise model.

% IVM

model.c(index, :) = 1./sqrt(model.varSigma(index, :));
for i = index
  model.u(i, :) = model.c(i, :).*(model.mu(i, :) ...
				  + model.noise.bias - model.y(i, :));
end

etaPart = model.noise.eta/(1-2*model.noise.eta);

for c = 1:size(model.y, 2)
  
  % Updates for the first category.
  indexStart = find(model.y(index, c) == 0);
  if ~isempty(indexStart)
    
    halfc = model.c(indexStart, c)/2;
    uMinusHalfc =  model.u(indexStart, c)-halfc;
    model.g(indexStart, c) = -model.c(indexStart, c)...
	.*ngaussian(uMinusHalfc) ...
	./(1-cummGaussian(uMinusHalfc) ...
	   + etaPart);
    
    model.nu(indexStart, c) = model.g(indexStart, c).*(model.g(indexStart, c) + model.c(indexStart, c).*uMinusHalfc);
  end
  
  % Updates for the last category.
  indexEnd = find(model.y(index, c) == model.C(c)-1);
  if ~isempty(indexEnd)
    
    halfc = model.c(indexEnd, c)/2;
    uPlusHalfc = model.u(indexEnd, c)+halfc;
    model.g(indexEnd, c) = model.c(indexEnd, c)...
	.*ngaussian(uPlusHalfc) ...
	./(cummGaussian(uPlusHalfc) ...
	   + etaPart);
    model.nu(indexEnd, c) = model.g(indexEnd, c).*(model.g(indexEnd, c) + model.c(indexEnd, c).*uPlusHalfc);
  end
  
  % Updates for other categories
  indexMiddle = find(model.y(index, c) ~= 0 & model.y(index, c) ~= model.C(c)-1);
  if ~isempty(indexMiddle)
    
    halfc = model.c(indexMiddle, c)/2;
    uPlusHalfc = model.u(indexMiddle, c)+halfc;
    uMinusHalfc =  model.u(indexMiddle, c)-halfc;
    model.g(indexMiddle, c) = model.c(indexMiddle, c)...
	.*(ngaussian(uPlusHalfc) ...
	   - ngaussian(uMinusHalfc)) ...
	./(cummGaussian(uPlusHalfc) ...
	   - cummGaussian(uMinusHalfc) ...
	   + etaPart);  
    
    model.nu(indexMiddle, c) = model.g(indexMiddle, c).*model.g(indexMiddle, c) ...
	+ model.c(indexMiddle, c).*model.c(indexMiddle, c)...
	.*(uPlusHalfc.*ngaussian(uPlusHalfc) ...
	   - uMinusHalfc.*ngaussian(uMinusHalfc))...
	./(cummGaussian(uPlusHalfc) ...
	   - cummGaussian(uMinusHalfc) ...
	   + etaPart);
  end
  model.nu(index, c) = thetaConstrain(model.nu(index, c));
end
