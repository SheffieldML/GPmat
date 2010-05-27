function K = sdlfmvXsdlfmKernComputeBlock(lfmKern1, lfmKern2, t1, t2, ...
    kyy, kyv, kvy, kvv, i, j, generalConst)

% SDLFMVXSDLFMKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% FORMAT
% DESC computes the kernel matrix for the SDLFM kernel function in the
% block specified at indeces i,j. It assumes the computation for functions
% that describe velocity (system 1) and position (system 2)
% ARG lfmKern1 : structure containing parameters for the system 1
% ARG lfmKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% ARG kyv : covariance for the initial conditions between position 1 and
% velocity 2 at block i,j
% ARG kvy : covariance for the initial conditions between velocity 1 and
% position 2 at block i,j
% ARG kvv : covariance for the initial conditions between velocity 1 and
% velocity 2 at block i,j
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdlfmKernComputeConstant.m
% RETURN K : the kernel matrix portion of block i,j
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin<11
    j = i;
    generalConst = [];
end

g1 = sdlfmvMeanCompute(lfmKern1(1), t1, 'Pos');
h1 = sdlfmvMeanCompute(lfmKern1(1), t1, 'Vel');
c2 = sdlfmMeanCompute(lfmKern2(1), t2, 'Pos');
e2 = sdlfmMeanCompute(lfmKern2(1), t2, 'Vel');

K = kyy*g1*c2.' + kyv*g1*e2.' + kvy*h1*c2.' + kvv*h1*e2.';

if i==j
    for k=1:length(lfmKern1)
        K  = K + lfmvXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, t2);
    end
else    
    if i>j
        PosPos = zeros(1, length(t2));
        PosVel = zeros(1, length(t2));
        for k =1:length(lfmKern1)
            PosPos = PosPos + lfmXlfmKernCompute(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
        end
        if isempty(generalConst{i,j})
            K = K + g1*PosPos + h1*PosVel;
        else
           
            K = K + (generalConst{i,j}(1,1)*g1 + generalConst{i,j}(2,1)*h1)*PosPos + ...
                (generalConst{i,j}(1,2)*g1 + generalConst{i,j}(2,2)*h1)*PosVel;
        end
    else
        PosVel = zeros(length(t1), 1);
        VelVel = zeros(length(t1), 1);
        for k =1:length(lfmKern1)
            PosVel = PosVel + lfmvXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
            VelVel = VelVel + lfmvXlfmvKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
        end
        if isempty(generalConst{i,j})
            K = K + PosVel*c2.' + VelVel*e2.';
        else
            K = K + PosVel*(generalConst{i,j}(1,1)*c2.' + generalConst{i,j}(2,1)*e2.') + ...
                VelVel*(generalConst{i,j}(1,2)*c2.' + generalConst{i,j}(2,2)*e2.');
        end
    end
end

