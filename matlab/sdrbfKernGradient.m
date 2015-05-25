function g = sdrbfKernGradient(sdrbfKern, t1, t2, covGrad)

% SDRBFKERNGRADIENT Gradient of SDRBF kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the switching
% dynamical radial basis function kernel's parameters. As well as the 
% kernel structure and the input positions, the user provides a matrix 
% PARTIAL which gives the partial derivatives of the function with respect 
% to the relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG t : the input locations for which the gradients are being
% computed. 
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in t1.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG t1 : the input locations associated with the rows of the
% kernel matrix.
% ARG t2 : the input locations associated with the columns of the
% kernel matrix.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as t1 and the same number of columns
% as t2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO sdrbfKernParamInit, kernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    covGrad = t2;
    t2 = t1;
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

%Form the basic kernels
rbfKern = struct();

% Create structures that will make easy the computation of the kernels
% spVector = [sdrbfKern.switchingTimes t1(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdrbfKern.switchingTimes) t1(end)+100];

dim1 = zeros(1, sdrbfKern.nIntervals); 
dim2 = zeros(1, sdrbfKern.nIntervals);

for i =1:sdrbfKern.nIntervals
    for j=1:sdrbfKern.nlfPerInt
        rbfKern(i,j).variance = sdrbfKern.variance;
        rbfKern(i,j).inverseWidth = sdrbfKern.inverseWidth(j,i);
        rbfKern(i,j).isNormalised = sdrbfKern.isNormalised;
    end
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2 = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

gradIW = zeros(sdrbfKern.nlfPerInt, sdrbfKern.nIntervals);
gradSP = zeros(1, sdrbfKern.nIntervals);
%dimParam = sdrbfKern.nIntervals*ones(1, 2);
%gC = cell(1, sdrbfKern.nlfPerInt);
%indxIW = 1:dimParam(1);
%indxSP = dimParam(1)+1:sum(dimParam);

covGrad2 = covGrad;
for q=1:sdrbfKern.nlfPerInt
    startValOne = 1;
    endValOne   = 0;
    startValThree = 1;
    endValThree = 0;
    gradSPTemp = zeros(1, sdrbfKern.nIntervals);
    if iscell(covGrad2)
        covGrad = covGrad2{q};
    end
    for i=1:sdrbfKern.nIntervals
        endValOne = endValOne + dim1(i);
        endValThree = endValThree + dim2(i);
        gP = rbfKernGradient(rbfKern(i,q), t1(startValOne:endValOne) - spVector(i), ...
            t2(startValThree:endValThree) - spVector(i), ...
            covGrad(startValOne:endValOne, startValThree:endValThree));
        gradIW(q,i)= gP(1);
        gX = rbfKernGradX(rbfKern(i,q), t1(startValOne:endValOne) - spVector(i), ...
            t2(startValThree:endValThree) - spVector(i));
        gXR = reshape(gX, dim1(i), dim2(i))';
        gradSPTemp(i) = gradSPTemp(i) - 2*sum(sum(gXR.*covGrad(startValOne:endValOne, startValThree:endValThree)));
        startValThree = endValThree + 1;
        startValOne = endValOne + 1;
    end
    gradSP = gradSP + gradSPTemp;
    %gC{q}(indxIW) = gradIW(q,:);
    %gC{q}(indxSP) = gradSPTemp; 
end

%g = [gradIW(:)' gradSP]; 

%%%%%%
tempGradSP = fliplr(gradSP);
tempGradSP = cumsum(tempGradSP); 
g = [gradIW(:)' fliplr(tempGradSP)]; 
%%%%%%
