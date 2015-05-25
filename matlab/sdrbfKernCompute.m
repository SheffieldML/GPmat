function [K, Kc] = sdrbfKernCompute(sdrbfKern, t1, t2)

% SDRBFKERNCOMPUTE Compute the SDRBF kernel given the parameters and t1.
% FORMAT
% DESC computes the kernel matrices for the switching dynamical radial
% basis function kernel. It is used together with the switching dynamical
% latent force model. The kernel matrices are stored in a cell structure
% and each component of the cell corresponds to the kernel matrix of a
% particular sequence of latent forces.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input points associated with the rows of the kernel.
% ARG t2 : the input points associated with the columns of the kernel.
% RETURN K : the kernel matrices computed at the given points. 
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations.
%
% FORMAT
% DESC computes the kernel matrices for the switching dynamical radial
% basis function kernel. It is used together with the switching dynamical
% latent force model. The kernel matrices are stored in a cell structure
% and each component of the cell corresponds to the kernel matrix of a
% particular sequence of latent forces.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : input data matrix in the form of a design matrix.
% RETURN K : the kernel matrix computed at the given points.
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations.
%
% SEEALSO : rbfKernParamInit, kernCompute, kernCreate, rbfKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN


if nargin < 3
    t2 = t1;
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

%Form the basic kernels
rbfKern = struct();

% Create structures that will make easy the computation of the kernels
% spVector = [sdrbfKern.switchingTimes t1(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdrbfKern.switchingTimes) t1(end)+50];

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

if sum(dim1)~=length(t1) || sum(dim2)~=length(t2)
    error('A problem with the dimensions of the switching intervals occured')
end

% Compute some necessary constants

K = cell(sdrbfKern.nlfPerInt,1);
Kc = cell(sdrbfKern.nlfPerInt,1);


for q=1:sdrbfKern.nlfPerInt
    startValOne = 1;
    endValOne   = 0;
    startValThree = 1;
    endValThree = 0;        
    tK = zeros(sum(dim1), sum(dim2));
    for i=1:sdrbfKern.nIntervals
        endValOne = endValOne + dim1(i);
        endValThree = endValThree + dim2(i);
        tK(startValOne:endValOne, startValThree:endValThree) = rbfKernCompute(rbfKern(i,q),...
            t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(i));
        Kc{q}{i} = tK(startValOne:endValOne, startValThree:endValThree);
        startValThree = endValThree + 1;
        startValOne = endValOne + 1;
    end
    tK = real(tK);
    K{q} = tK;
end

