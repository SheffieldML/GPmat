function [g1, g2, g3] = sdlfmXsdlfmKernGradientBlockIGJ(lfmKern1, lfmKern2, ...
    t1, t2,  i, j, generalConst, generalConstGrad, covGrad, g1Mean, ...
    g2Mean, gsp1Mean, gsp2Mean, typeMeanParam, typeMeanSwitching, ...
    typeKernParam1, typeKernParam2, typeKernSwitching1, typeKernSwitching2)

% SDLFMXSDLFMKERNGRADIENTBLOCKIGJ 
% FORMAT
% DESC computes the gradients of the parameters for system 1 and system 2
% when i is greater than j.
% ARG lfmKern1 : structure containing parameters for the system 1
% ARG lfmKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdlfmKernComputeConstant.m
% ARG generalConstGrad : derivatives of the constants computed with
% sdlfmKernGradientConstant.m
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the corresponding kernel matrix
% ARG g1Mean : gradients of the parameter of the system 1 obtained from the
% part of the function that uses the funcitons accommpanying the initial
% conditions.
% ARG g2Mean : gradients of the parameter of the system 2 obtained from the
% part of the function that uses the funcitons accommpanying the initial
% conditions.
% ARG gsp1Mean : gradient of the switching point of the system 1 obtained 
% from the part of the function that uses the funcitons accommpanying the 
% initial conditions.
% ARG gsp2Mean : gradient of the switching point of the system 2 obtained 
% from the part of the function that uses the funcitons accommpanying the 
% initial conditions.
% ARG typeMeanParam : specify the mean functions used to compute this part 
% of the kernel
% ARG typeMeanSwitching : specify the functions used to compute the
% gradients of the swicthing points in part thst uses mean functions
% ARG typeKernParam1 : specify the first kernel function used to compute 
% this part of the kernel
% ARG typeKernParam2 : specify the second kernel function used to compute 
% this part of the kernel
% ARG typeKernSwitching1 : specify the functions used to compute the
% gradients of the swicthing points in both sides of the kernel function 1
% ARG typeKernSwitching2 : specify the functions used to compute the
% gradients of the swicthing points in both sides of the kernel function 2
% RETURN g1 : gradients of parameters for the system 1
% RETURN g2 : gradients of parameters for the system 2
% RETURN g3 : gradients of switching points
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin < 14
    typeMeanParam = 'sdlfm';
    typeMeanSwitching = 'sdlfmv';
    typeKernParam1 = 'lfmXlfm';
    typeKernParam2 = 'lfmvXlfm';
    typeKernSwitching1= {'lfmvXlfm', 'lfmvXlfm'};
    typeKernSwitching2 ={'lfmaXlfm', 'lfmvXlfmv'};
end

fhandleMeanParam = str2func([typeMeanParam 'MeanCompute']);
fhandleMeanGradParam = str2func([typeMeanParam 'MeanGradient']);
fhandleMeanGradSwitching = str2func([typeMeanSwitching 'MeanCompute']);
fhandleKernPosPos = str2func([typeKernParam1 'KernCompute']);
fhandleKernGradPosPos = str2func([typeKernParam1 'KernGradient']);
fhandleKernGradSwitchingPosPos1 = str2func([typeKernSwitching1{1} 'KernCompute']);
fhandleKernGradSwitchingPosPos2 = str2func([typeKernSwitching1{2} 'KernCompute']);
fhandleKernVelPos = str2func([typeKernParam2 'KernCompute']);
fhandleKernGradVelPos = str2func([typeKernParam2 'KernGradient']);
fhandleKernGradSwitchingVelPos1 = str2func([typeKernSwitching2{1} 'KernCompute']);
fhandleKernGradSwitchingVelPos2 = str2func([typeKernSwitching2{2} 'KernCompute']);


c1 = fhandleMeanParam(lfmKern1(1), t1, 'Pos');
e1 = fhandleMeanParam(lfmKern1(1), t1, 'Vel');

if isempty(generalConst{i,j})
    coeffPosPos = c1;
    coeffVelPos = e1;
else
    coeffPosPos = generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1;
    coeffVelPos = generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1;
end
g1Kern = zeros(length(lfmKern1), 5);
g2Kern = zeros(length(lfmKern1), 5);
PosPos = zeros(1,length(t2));
VelPos = zeros(1,length(t2));

for k=1:length(lfmKern1)
    PosPos = PosPos + fhandleKernPosPos(lfmKern1(k), lfmKern2(k), ...
        lfmKern2(k).limit, t2);
    [g1KLocal, g2KLocal] = fhandleKernGradPosPos(lfmKern1(k), lfmKern2(k), ...
        lfmKern2(k).limit, t2, covGrad.', coeffPosPos);
    g1Kern(k,:) = g1Kern(k, :) + g1KLocal;
    g2Kern(k,:) = g2Kern(k, :) + g2KLocal;
    VelPos = VelPos + fhandleKernVelPos(lfmKern1(k), lfmKern2(k), ...
        lfmKern2(k).limit, t2);
    [g1KLocal, g2KLocal] = fhandleKernGradVelPos(lfmKern1(k), lfmKern2(k), ...
        lfmKern2(k).limit, t2, covGrad.', coeffVelPos);
    g1Kern(k,:) = g1Kern(k, :) + g1KLocal;
    g2Kern(k,:) = g2Kern(k, :) + g2KLocal;
end
[gcAlpha, gcOmega] = fhandleMeanGradParam(lfmKern1(1), t1, 'Pos');
[geAlpha, geOmega] = fhandleMeanGradParam(lfmKern1(1), t1, 'Vel');
[gradAlpha, gradOmega] = getLocalGradAlphaOmega(lfmKern1);
g1Pos = fhandleMeanGradSwitching(lfmKern1(1), t1, 'Pos');
h1Vel = fhandleMeanGradSwitching(lfmKern1(1), t1, 'Vel');
gsp1Kern = zeros(1, length(lfmKern1));
gsp2Kern = zeros(1, length(lfmKern1));
% This means that are only two switching points involved, t1
% and t0 appear in k, like k_{ff}(t1-t0,t'-t0) and
% k_{mf}(t1-t0,t'-t0). The other points appear in c1(t - t1)
% and e1(t - t1)
for k=1:length(lfmKern1)
    temp = fhandleKernGradSwitchingPosPos1(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
    gsp2Kern(k) = gsp2Kern(k) - sum(sum((coeffPosPos*temp).*covGrad));
    gsp1Kern(k) = gsp1Kern(k) + sum(sum((coeffPosPos*temp).*covGrad));
    temp = fhandleKernGradSwitchingPosPos2(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).';
    gsp2Kern(k) = gsp2Kern(k) - sum(sum((coeffPosPos*temp).*covGrad));
    temp = fhandleKernGradSwitchingVelPos1(lfmKern1(k), lfmKern2(k), lfmKern2(k).limit, t2);
    gsp2Kern(k) = gsp2Kern(k) - sum(sum((coeffVelPos*temp).*covGrad));
    gsp1Kern(k) = gsp1Kern(k) + sum(sum((coeffVelPos*temp).*covGrad));
    temp = fhandleKernGradSwitchingVelPos2(lfmKern2(k), lfmKern1(k), t2, lfmKern2(k).limit).';
    gsp2Kern(k) = gsp2Kern(k) - sum(sum((coeffVelPos*temp).*covGrad));
end
if isempty(generalConst{i,j})
    matGradAlpha = gcAlpha*PosPos + geAlpha*VelPos;
    matGradOmega = gcOmega*PosPos + geOmega*VelPos;
    gCoeff = gradAlpha*sum(sum(matGradAlpha.*covGrad)) + ...
        gradOmega*(sum(sum(matGradOmega.*covGrad)));
    matGradSp1 = g1Pos*PosPos + h1Vel*VelPos;
    gsp1 = - sum(sum(matGradSp1.*covGrad));
    g3(i) = gsp1Mean + sum(gsp1Kern) + gsp1; % switching point 1
    g3(j) = gsp2Mean + sum(gsp2Kern);        % switching point 2
else
    constGradAlpha = generalConstGrad{1}{i,j};
    constGradOmega = generalConstGrad{2}{i,j};
    constVal = generalConst{i,j};
    % Derivative wrt parameters
    matGradAlpha = (constVal(1,1)*gcAlpha + constGradAlpha(1,1)*c1 ...
        + constVal(2,1)*geAlpha + constGradAlpha(2,1)*e1)*PosPos ...
        + (constVal(1,2)*gcAlpha + constGradAlpha(1,2)*c1 ...
        + constVal(2,2)*geAlpha + constGradAlpha(2,2)*e1)*VelPos;
    matGradOmega = (constVal(1,1)*gcOmega + constGradOmega(1,1)*c1 ...
        + constVal(2,1)*geOmega + constGradOmega(2,1)*e1)*PosPos ...
        + (constVal(1,2)*gcOmega + constGradOmega(1,2)*c1 ...
        + constVal(2,2)*geOmega + constGradOmega(2,2)*e1)*VelPos;
    gCoeff = gradAlpha*sum(sum(matGradAlpha.*covGrad)) + ...
        gradOmega*(sum(sum(matGradOmega.*covGrad)));
    % Derivative wrt switching points
    % Firts, notice that gsp1Kern doesn't correspond to the switching point
    % of the interval i (or say 1). It's just the derivative of the inner
    % switching point that lead to the actual switching point for the
    % interval 1. So we rename it first.
    gspInt = sum(gsp1Kern);
    % Compute the derivative of the swicthing point i, not appearing in the
    % constant
    matGradSp1 = (constVal(1,1)*g1Pos + constVal(2,1)*h1Vel)*PosPos + ...
        (constVal(1,2)*g1Pos + constVal(2,2)*h1Vel)*VelPos;
    gsp1 = - sum(sum(matGradSp1.*covGrad));
    % Compute the derivatives for all the switching points appearing in the
    % constant. The order is descending, i.e., t_3_, t_2, t_1
    constGradSPoint = generalConstGrad{3}{i,j};
    numberSP = size(constGradSPoint,2);
    gspInBetween = zeros(1, numberSP);
    for k=1:numberSP
        temp = (constGradSPoint(1,k)*c1 + constGradSPoint(3,k)*e1)*PosPos + ...
            (constGradSPoint(2,k)*c1 + constGradSPoint(4,k)*e1)*VelPos;
        gspInBetween(k) = sum(sum(temp.*covGrad));
    end
    % Assign derivatives wrt all other switching points, with a correction
    % for the derivative of the innermost switching point in the constant
    gspInBetween(end) =  gspInBetween(end) + gspInt;
    g3(i) = gsp1Mean + gsp1 + gspInBetween(1);
    g3(j) = gsp2Mean + sum(gsp2Kern);
    g3(j+1:i-1) = fliplr(gspInBetween(2:end));
end
% Assign derivatives wrt first system
g1{1} = g1Mean + sum(g1Kern(:,1:3), 1) + gCoeff; % mass 1, spring 1, damper 1
g1{2} = g1Kern(:,4).';                            % inverse widths
g1{3} = g1Kern(:,5).';                            % seinsitivities 1
% Assign derivatives wrt second system
g2{1} = g2Mean + sum(g2Kern(:,1:3), 1); % mass 2, spring 2, damper 2
g2{2} = g2Kern(:,4).';                   % inverse widths
g2{3} = g2Kern(:,5).';


