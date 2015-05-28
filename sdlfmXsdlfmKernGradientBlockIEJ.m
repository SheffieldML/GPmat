function [g1, g2, g3] = sdlfmXsdlfmKernGradientBlockIEJ(lfmKern1, ...
    lfmKern2, t1, t2, covGrad,  g1Mean, g2Mean, gsp1Mean, gsp2Mean, ...
    typeParam, typeSwitching1, typeSwitching2)

% SDLFMXSDLFMKERNGRADIENTBLOCKIEJ 
% FORMAT
% DESC computes the gradients of the parameters for system 1 and system 2
% when i is equal to j, this is, when the kernel function is evaluated
% at the same switching interval
% ARG lfmKern1 : structure containing parameters for the system 1
% ARG lfmKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
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
% ARG typeParam : specify the mean functions used to compute this part of
% the kernel
% ARG typeSwitching1 : specify the functions used to compute the
% gradients of the swicthing points in the left side of the kernel
% ARG typeSwitching2 : specify the functions used to compute the
% gradients of the swicthing points in the right side of the kernel
% RETURN g1 : gradients of parameters for the system 1
% RETURN g2 : gradients of parameters for the system 2
% RETURN g3 : gradients of switching points
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin < 10
    typeParam{1} = 'lfm';
    typeParam{2} = 'lfm';
    typeSwitching1{1} = 'lfmv';
    typeSwitching1{2} = 'lfm';
    typeSwitching2{1} = 'lfmv';
    typeSwitching2{2} = 'lfm';
end

fhandleGradParam = str2func([typeParam{1} 'X' typeParam{2} 'KernGradient']);
fhandleGradSPoint1 = str2func([typeSwitching1{1} 'X' typeSwitching1{2} 'KernCompute']);    
fhandleGradSPoint2 = str2func([typeSwitching2{1} 'X' typeSwitching2{2} 'KernCompute']);    

g1Kern = zeros(length(lfmKern1), 5);
g2Kern = zeros(length(lfmKern1), 5);
gsp1Kern = zeros(1, length(lfmKern1));
gsp2Kern = zeros(1, length(lfmKern1));
for k=1:length(lfmKern1)
    % Compute the derivatives of the kernel function wrt parameters
    [g1Kern(k, :), g2Kern(k, :)] = fhandleGradParam(lfmKern1(k), ...
        lfmKern2(k), t1, t2, covGrad);
    % Compute derivative of the switching points
    temp = fhandleGradSPoint1(lfmKern1(k), lfmKern2(k), t1, t2);
    gsp1Kern(k) = - sum(sum(temp.*covGrad));
    temp = fhandleGradSPoint2(lfmKern2(k), lfmKern1(k), t2, t1).';
    gsp2Kern(k) = - sum(sum(temp.*covGrad));
end
% Assign derivatives wrt first system
g1{1} = g1Mean + sum(g1Kern(:,1:3), 1); % mass 1, spring 1, damper 1
g1{2} = g1Kern(:,4).';                   % inverse widths
g1{3} = g1Kern(:,5).';                   % sensitivities 1
% Assign derivatives wrt second system
g2{1} = g2Mean + sum(g2Kern(:,1:3), 1); % mass 2, spring 2, damper 2
g2{2} = g2Kern(:,4).';                   % inverse widths
g2{3} = g2Kern(:,5).';                   % sensitivities 2

g3 = gsp1Mean + sum(gsp1Kern) + gsp2Mean + sum(gsp2Kern); % All these are the same switching point.

