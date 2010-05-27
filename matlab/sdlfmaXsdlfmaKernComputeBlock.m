function K = sdlfmaXsdlfmaKernComputeBlock(lfmKern1, lfmKern2, t1, t2, ...
    kyy, kyv, kvy, kvv, i, j, generalConst)

% SDLFMAXSDLFMAKERNCOMPUTEBLOCK Computes SDLFM kernel matrix for block i,j
% FORMAT
% DESC computes the kernel matrix for the SDLFM kernel function in the
% block specified at indeces i,j. It assumes the computation for functions
% that describe accelerations (acceleration 1 and acceleration 2).
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

a1 = sdlfmaMeanCompute(lfmKern1(1), t1, 'Pos');
b1 = sdlfmaMeanCompute(lfmKern1(1), t1, 'Vel');
a2 = sdlfmaMeanCompute(lfmKern2(1), t2, 'Pos');
b2 = sdlfmaMeanCompute(lfmKern2(1), t2, 'Vel');

K = kyy*a1*a2.' + kyv*a1*b2.' + kvy*b1*a2.' + kvv*b1*b2.';

if i==j
    for k=1:length(lfmKern1)
        K  = K + lfmaXlfmaKernCompute(lfmKern1(k), lfmKern2(k), t1, t2);
    end
else    
    if i>j
        AccelPos = zeros(1, length(t2));
        AccelVel = zeros(1, length(t2));
        for k=1:length(lfmKern1)
            AccelPos = AccelPos + lfmaXlfmKernCompute(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).'; 
            AccelVel = AccelVel + lfmaXlfmvKernCompute(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).';
        end
        if isempty(generalConst{i,j})
            K = K + a1*AccelPos + b1*AccelVel;        
        else
            K = K + (generalConst{i,j}(1,1)*a1 + generalConst{i,j}(2,1)*b1)*AccelPos + ...
                (generalConst{i,j}(1,2)*a1 + generalConst{i,j}(2,2)*b1)*AccelVel;           
        end 
    else
        AccelPos = zeros(length(t1),1);
        AccelVel = zeros(length(t1),1);
        for k =1:length(lfmKern1)
            AccelPos = AccelPos + lfmaXlfmKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
            AccelVel = AccelVel + lfmaXlfmvKernCompute(lfmKern1(k), lfmKern2(k), t1, lfmKern1(k).limit);
        end
        if isempty(generalConst{i,j})
            K = K + AccelPos*a2.' + AccelVel*b2.';
        else
            K = K + AccelPos*(generalConst{i,j}(1,1)*a2.' + generalConst{i,j}(2,1)*b2.') + ...
                AccelVel*(generalConst{i,j}(1,2)*a2.' + generalConst{i,j}(2,2)*b2.');
        end
    end
end

