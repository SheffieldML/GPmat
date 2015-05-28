function [generalConstGrad, generalConst] = sdlfmKernGradientConstant(nIntervals, ...
    lfmKern1, lfmKern2, spVector)

% SDLFMKERNGRADIENTCONSTANT Gradients for constants for the SDLFM kernel
% FORMAT
% DESC computes necessary gradients of the constants computed with
% sdlfmKernComputeConstant.m
% ARG nIntervals : number of switching intervals in the kernel
% ARG lfmKern1 : structure containing the parameters of system 1
% ARG lfmkern2 : structure containing the parameters of system 2
% ARG spVector : vector containing the switching time values
% RETURN generalConstGrad : a cell containing the relative derivatives
% RETURN generalConstant : a cell containing the necessary constants for
% computing the kernel in the switching intervals.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

generalConst = cell(nIntervals);
generalGradAlpha = cell(nIntervals);
generalGradOmega = cell(nIntervals);
generalGradSPoint = cell(nIntervals);

for i=1:nIntervals
    for j =1:i-1
        if i - j>=2
            c1 = sdlfmMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            e1 = sdlfmMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            c2 = sdlfmMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            e2 = sdlfmMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            g1 = sdlfmvMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            h1 = sdlfmvMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            g2 = sdlfmvMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            h2 = sdlfmvMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            % Derivatives wrt alpha and omega so that we can get the
            % derivatives wrt paramters in a further step.
            [gc1Alpha, gc1Omega] = sdlfmMeanGradient(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            [ge1Alpha, ge1Omega] = sdlfmMeanGradient(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            [gc2Alpha, gc2Omega] = sdlfmMeanGradient(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            [ge2Alpha, ge2Omega] = sdlfmMeanGradient(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            [gg1Alpha, gg1Omega] = sdlfmvMeanGradient(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            [gh1Alpha, gh1Omega] = sdlfmvMeanGradient(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            [gg2Alpha, gg2Omega] = sdlfmvMeanGradient(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            [gh2Alpha, gh2Omega] = sdlfmvMeanGradient(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            % Derivatives wrt the switching points
            gc1sp = sdlfmvMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            ge1sp = sdlfmvMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            gc2sp = sdlfmvMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            ge2sp = sdlfmvMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            gg1sp = sdlfmaMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Pos');
            gh1sp = sdlfmaMeanCompute(lfmKern1, spVector(i) - spVector(i-1), 'Vel');
            gg2sp = sdlfmaMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Pos');
            gh2sp = sdlfmaMeanCompute(lfmKern2, spVector(i) - spVector(i-1), 'Vel');
            if i - j == 2
                % Computation of the constants
                generalConst{i,j}(1,1) = c1; generalConst{i,j}(1,2) = e1;
                generalConst{i,j}(2,1) = g1; generalConst{i,j}(2,2) = h1;
                generalConst{j,i}(1,1) = c2; generalConst{j,i}(1,2) = e2;
                generalConst{j,i}(2,1) = g2; generalConst{j,i}(2,2) = h2;
                % Gradient wrt alpha for lfmKern1
                generalGradAlpha{i,j}(1,1) = gc1Alpha; generalGradAlpha{i,j}(1,2) = ge1Alpha;
                generalGradAlpha{i,j}(2,1) = gg1Alpha; generalGradAlpha{i,j}(2,2) = gh1Alpha;
                % Gradient wrt alpha for lfmKern2
                generalGradAlpha{j,i}(1,1) = gc2Alpha; generalGradAlpha{j,i}(1,2) = ge2Alpha;
                generalGradAlpha{j,i}(2,1) = gg2Alpha; generalGradAlpha{j,i}(2,2) = gh2Alpha;
                % Gradient wrt omega for lfmKern1
                generalGradOmega{i,j}(1,1) = gc1Omega; generalGradOmega{i,j}(1,2) = ge1Omega;
                generalGradOmega{i,j}(2,1) = gg1Omega; generalGradOmega{i,j}(2,2) = gh1Omega;
                % Gradient wrt omega for lfmKern2
                generalGradOmega{j,i}(1,1) = gc2Omega; generalGradOmega{j,i}(1,2) = ge2Omega;
                generalGradOmega{j,i}(2,1) = gg2Omega; generalGradOmega{j,i}(2,2) = gh2Omega;
                % Gradient wrt switching points for lfmKern1
                generalGradSPoint{i,j}(1,1) = gc1sp; generalGradSPoint{i,j}(1,2) = -gc1sp;
                generalGradSPoint{i,j}(2,1) = ge1sp; generalGradSPoint{i,j}(2,2) = -ge1sp;
                generalGradSPoint{i,j}(3,1) = gg1sp; generalGradSPoint{i,j}(3,2) = -gg1sp;
                generalGradSPoint{i,j}(4,1) = gh1sp; generalGradSPoint{i,j}(4,2) = -gh1sp;
                % Gradient wrt switching points for lfmKern2
                generalGradSPoint{j,i}(1,1) = gc2sp; generalGradSPoint{j,i}(1,2) = -gc2sp; % f1 wrt t_i and t_{i-1}
                generalGradSPoint{j,i}(2,1) = ge2sp; generalGradSPoint{j,i}(2,2) = -ge2sp; % f2 wrt t_i and t_{i-1}
                generalGradSPoint{j,i}(3,1) = gg2sp; generalGradSPoint{j,i}(3,2) = -gg2sp; % f3 wrt t_i and t_{i-1}
                generalGradSPoint{j,i}(4,1) = gh2sp; generalGradSPoint{j,i}(4,2) = -gh2sp; % f4 wrt t_i and t_{i-1}
                
            else
                % Derivatives wrt to lfmKern1
                % Allocate values to make callings to
                % computeLocalDerivative short
                f1 = generalConst{i-1,j}(1,1); f2 = generalConst{i-1,j}(1,2);
                f3 = generalConst{i-1,j}(2,1); f4 = generalConst{i-1,j}(2,2);
                gf1Alpha = generalGradAlpha{i-1,j}(1,1); gf2Alpha = generalGradAlpha{i-1,j}(1,2);
                gf3Alpha = generalGradAlpha{i-1,j}(2,1); gf4Alpha = generalGradAlpha{i-1,j}(2,2);
                gf1Omega = generalGradOmega{i-1,j}(1,1); gf2Omega = generalGradOmega{i-1,j}(1,2);
                gf3Omega = generalGradOmega{i-1,j}(2,1); gf4Omega = generalGradOmega{i-1,j}(2,2);                
                % Compute constants related to lfmKern1
                generalConst{i,j}(1,1) = c1*generalConst{i-1,j}(1,1) + e1*generalConst{i-1,j}(2,1);
                generalConst{i,j}(1,2) = c1*generalConst{i-1,j}(1,2) + e1*generalConst{i-1,j}(2,2);
                generalConst{i,j}(2,1) = g1*generalConst{i-1,j}(1,1) + h1*generalConst{i-1,j}(2,1);
                generalConst{i,j}(2,2) = g1*generalConst{i-1,j}(1,2) + h1*generalConst{i-1,j}(2,2);
                % Gradient wrt alpha for lfmKern1
                generalGradAlpha{i,j}(1,1) = computeLocalDerivative(c1, f1, e1, f3, gc1Alpha, gf1Alpha, ge1Alpha, gf3Alpha);
                generalGradAlpha{i,j}(1,2) = computeLocalDerivative(c1, f2, e1, f4, gc1Alpha, gf2Alpha, ge1Alpha, gf4Alpha);
                generalGradAlpha{i,j}(2,1) = computeLocalDerivative(g1, f1, h1, f3, gg1Alpha, gf1Alpha, gh1Alpha, gf3Alpha);
                generalGradAlpha{i,j}(2,2) = computeLocalDerivative(g1, f2, h1, f4, gg1Alpha, gf2Alpha, gh1Alpha, gf4Alpha);
                % Gradient wrt omega for lfmKern1
                generalGradOmega{i,j}(1,1) = computeLocalDerivative(c1, f1, e1, f3, gc1Omega, gf1Omega, ge1Omega, gf3Omega);
                generalGradOmega{i,j}(1,2) = computeLocalDerivative(c1, f2, e1, f4, gc1Omega, gf2Omega, ge1Omega, gf4Omega);
                generalGradOmega{i,j}(2,1) = computeLocalDerivative(g1, f1, h1, f3, gg1Omega, gf1Omega, gh1Omega, gf3Omega);
                generalGradOmega{i,j}(2,2) = computeLocalDerivative(g1, f2, h1, f4, gg1Omega, gf2Omega, gh1Omega, gf4Omega);
                % Gradient wrt to switching point t_i and t_{i-1}
                gf1sp = generalGradSPoint{i-1,j}(1,1); gf2sp = generalGradSPoint{i-1,j}(2,1);
                gf3sp = generalGradSPoint{i-1,j}(3,1); gf4sp = generalGradSPoint{i-1,j}(4,1);
                % First the derivative wrt to t_i
                generalGradSPoint{i,j}(1,1) = gc1sp*f1 + ge1sp*f3;
                generalGradSPoint{i,j}(2,1) = gc1sp*f2 + ge1sp*f4;
                generalGradSPoint{i,j}(3,1) = gg1sp*f1 + gh1sp*f3;
                generalGradSPoint{i,j}(4,1) = gg1sp*f2 + gh1sp*f4;
                % Second the derivative wrt to t_{i-1}
                generalGradSPoint{i,j}(1,2) = computeLocalDerivative(c1, f1, e1, f3, -gc1sp, gf1sp, -ge1sp, gf3sp);
                generalGradSPoint{i,j}(2,2) = computeLocalDerivative(c1, f2, e1, f4, -gc1sp, gf2sp, -ge1sp, gf4sp);
                generalGradSPoint{i,j}(3,2) = computeLocalDerivative(g1, f1, h1, f3, -gg1sp, gf1sp, -gh1sp, gf3sp);
                generalGradSPoint{i,j}(4,2) = computeLocalDerivative(g1, f2, h1, f4, -gg1sp, gf2sp, -gh1sp, gf4sp);
                % Now compute all other derivatives wrt to the last
                % switching points
                maxL = i - j;
                cont = 1;
                for k = 3:maxL,
                    cont = cont + 1;
                    generalGradSPoint{i,j}(1,k) = c1*generalGradSPoint{i-1,j}(1,cont) + e1*generalGradSPoint{i-1,j}(3,cont);
                    generalGradSPoint{i,j}(2,k) = c1*generalGradSPoint{i-1,j}(2,cont) + e1*generalGradSPoint{i-1,j}(4,cont);
                    generalGradSPoint{i,j}(3,k) = g1*generalGradSPoint{i-1,j}(1,cont) + h1*generalGradSPoint{i-1,j}(3,cont);
                    generalGradSPoint{i,j}(4,k) = g1*generalGradSPoint{i-1,j}(2,cont) + h1*generalGradSPoint{i-1,j}(4,cont);
                end
                % Derivatives wrt to lfmKern2
                % Allocate values to make callings to
                % computeLocalDerivative short
                f1 = generalConst{j,i-1}(1,1); f2 = generalConst{j,i-1}(1,2);
                f3 = generalConst{j,i-1}(2,1); f4 = generalConst{j,i-1}(2,2);
                gf1Alpha = generalGradAlpha{j,i-1}(1,1); gf2Alpha = generalGradAlpha{j,i-1}(1,2);
                gf3Alpha = generalGradAlpha{j,i-1}(2,1); gf4Alpha = generalGradAlpha{j,i-1}(2,2);
                gf1Omega = generalGradOmega{j,i-1}(1,1); gf2Omega = generalGradOmega{j,i-1}(1,2);
                gf3Omega = generalGradOmega{j,i-1}(2,1); gf4Omega = generalGradOmega{j,i-1}(2,2);
                % Compute constants related to lfmKern2
                generalConst{j,i}(1,1) = c2*generalConst{j,i-1}(1,1) + e2*generalConst{j,i-1}(2,1);
                generalConst{j,i}(1,2) = c2*generalConst{j,i-1}(1,2) + e2*generalConst{j,i-1}(2,2);
                generalConst{j,i}(2,1) = g2*generalConst{j,i-1}(1,1) + h2*generalConst{j,i-1}(2,1);
                generalConst{j,i}(2,2) = g2*generalConst{j,i-1}(1,2) + h2*generalConst{j,i-1}(2,2);
                % Gradient wrt alpha for lfmKern2
                generalGradAlpha{j,i}(1,1) = computeLocalDerivative(c2, f1, e2, f3, gc2Alpha, gf1Alpha, ge2Alpha, gf3Alpha);
                generalGradAlpha{j,i}(1,2) = computeLocalDerivative(c2, f2, e2, f4, gc2Alpha, gf2Alpha, ge2Alpha, gf4Alpha);
                generalGradAlpha{j,i}(2,1) = computeLocalDerivative(g2, f1, h2, f3, gg2Alpha, gf1Alpha, gh2Alpha, gf3Alpha);
                generalGradAlpha{j,i}(2,2) = computeLocalDerivative(g2, f2, h2, f4, gg2Alpha, gf2Alpha, gh2Alpha, gf4Alpha);
                % Gradient wrt omega for lfmKern2
                generalGradOmega{j,i}(1,1) = computeLocalDerivative(c2, f1, e2, f3, gc2Omega, gf1Omega, ge2Omega, gf3Omega);
                generalGradOmega{j,i}(1,2) = computeLocalDerivative(c2, f2, e2, f4, gc2Omega, gf2Omega, ge2Omega, gf4Omega);
                generalGradOmega{j,i}(2,1) = computeLocalDerivative(g2, f1, h2, f3, gg2Omega, gf1Omega, gh2Omega, gf3Omega);
                generalGradOmega{j,i}(2,2) = computeLocalDerivative(g2, f2, h2, f4, gg2Omega, gf2Omega, gh2Omega, gf4Omega);
                % Gradient wrt to switching point t_i and t_{i-1}
                gf1sp = generalGradSPoint{j, i-1}(1,1); gf2sp = generalGradSPoint{j, i-1}(2,1);
                gf3sp = generalGradSPoint{j, i-1}(3,1); gf4sp = generalGradSPoint{j, i-1}(4,1);
                % First the derivative wrt to t_i
                generalGradSPoint{j,i}(1,1) = gc2sp*f1 + ge2sp*f3;
                generalGradSPoint{j,i}(2,1) = gc2sp*f2 + ge2sp*f4;
                generalGradSPoint{j,i}(3,1) = gg2sp*f1 + gh2sp*f3;
                generalGradSPoint{j,i}(4,1) = gg2sp*f2 + gh2sp*f4;
                % Second the derivative wrt to t_{i-1}
                generalGradSPoint{j,i}(1,2) = computeLocalDerivative(c2, f1, e2, f3, -gc2sp, gf1sp, -ge2sp, gf3sp);
                generalGradSPoint{j,i}(2,2) = computeLocalDerivative(c2, f2, e2, f4, -gc2sp, gf2sp, -ge2sp, gf4sp);
                generalGradSPoint{j,i}(3,2) = computeLocalDerivative(g2, f1, h2, f3, -gg2sp, gf1sp, -gh2sp, gf3sp);
                generalGradSPoint{j,i}(4,2) = computeLocalDerivative(g2, f2, h2, f4, -gg2sp, gf2sp, -gh2sp, gf4sp);
                % Now compute all other derivatives wrt to the last
                % switching points
                maxL = i - j;
                cont = 1;
                for k = 3:maxL,
                    cont = cont + 1;
                    generalGradSPoint{j,i}(1,k) = c2*generalGradSPoint{j,i-1}(1,cont) + e2*generalGradSPoint{j,i-1}(3,cont);
                    generalGradSPoint{j,i}(2,k) = c2*generalGradSPoint{j,i-1}(2,cont) + e2*generalGradSPoint{j,i-1}(4,cont);
                    generalGradSPoint{j,i}(3,k) = g2*generalGradSPoint{j,i-1}(1,cont) + h2*generalGradSPoint{j,i-1}(3,cont);
                    generalGradSPoint{j,i}(4,k) = g2*generalGradSPoint{j,i-1}(2,cont) + h2*generalGradSPoint{j,i-1}(4,cont);
                end
                
                
            end
        end
    end
end

generalConstGrad{1} = generalGradAlpha;
generalConstGrad{2} = generalGradOmega;
generalConstGrad{3} = generalGradSPoint;

function gradtot = computeLocalDerivative(a, b, c, d, grada, gradb, gradc, gradd)

gradtot = a*gradb + grada*b + c*gradd + gradc*d;
