function K = sdlfmXsdrbfKernComputeBlock(lfmKern, rbfKern, t1, t2, ...
    i, j, generalConst)

% SDLFMXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDLFM and SDRBF for i,j
% FORMAT
% DESC computes the kernel matrix between a SDLFM kernel function
% and a SDRBF kernel function in the block specified at indeces i,j. It 
% assumes the computation for a function that describe a position.
% ARG lfmKern : structure containing parameters for the outputs system
% ARG rbfKern : structure containing parameters for the latent system 
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdlfmKernComputeConstant.m
% RETURN K : the kernel matrix portion of block i,j
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin<7
    j = i;
    generalConst = [];
end

if i==j
    K  = lfmXrbfKernCompute(lfmKern, rbfKern, t1, t2);    
else
    c1 = sdlfmMeanCompute(lfmKern, t1, 'Pos');
    e1 = sdlfmMeanCompute(lfmKern, t1, 'Vel');
    if i>j
        PosRbf = lfmXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
        VelRbf = lfmvXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
        if isempty(generalConst{i,j})
            K =  c1*PosRbf + e1*VelRbf;
        else
            K = (generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1)*PosRbf + ...
                (generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1)*VelRbf;
        end
    end
end

