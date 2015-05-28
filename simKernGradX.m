function gT = simKernGradX(kern, t1, t2);

% SIMKERNGRADX Gradient of SIM kernel with respect to each time point in t1.
% FORMAT
% DESC computes the gradient of the single input motif
% kernel with respect to the input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG t : locations against which gradients are being computed.
% RETURN gT : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in t (currently always one).
%
% FORMAT
% DESC computes the gradident of the single input motif
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG t1 : row locations against which gradients are being computed.
% ARG t2 : column locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in t1, numData2 is the number of data
% points in t2 and numInputs is the number of input
% dimensions in t1 and t2 (currently always one).
%
% SEEALSO simKernParamInit, kernGradX, simKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : David Luengo, 2009

% KERN

if (size(t1, 2) > 1) || (size(t2, 2) > 1)
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);
sigma2 = sigma*sigma;
t1 = t1 - kern.delay;
t2 = t2 - kern.delay;
D = kern.decay;
halfSigmaD = 0.5*sigma*D;

gT = zeros(size(t2, 1), 1, size(t1, 1));

if (kern.isStationary == false)
    h1 = simComputeH(t1, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
    for i=1:size(t1, 1)
        diffT = t1(i)-t2;
        lnPart1 = lnDiffErfs(halfSigmaD+t1(i)/sigma, halfSigmaD+diffT/sigma);
        lnPart2 = lnDiffErfs(halfSigmaD, halfSigmaD-t2/sigma);
        gT(:, :, i) = 0.5 * kern.variance * (-D*h1(i,:).' ...
                + 1/(sigma*D*sqrt(pi))*(exp(-diffT.*diffT/sigma2)-exp(-D*t2-t1(i)*t1(i)/sigma2)) ...
            + 0.5 * (exp(halfSigmaD*halfSigmaD + D*diffT + lnPart1) ...
                + exp(halfSigmaD*halfSigmaD - D*(t1(i)+t2) + lnPart2)) ...
                + 1/(sigma*D*sqrt(pi)) * (exp(-D*t2-t1(i)*t1(i)/sigma2) - exp(-diffT.*diffT/sigma2)));
    end
else
    h1 = simComputeHStat(t1, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
    h2 = simComputeHStat(t2, t1, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
    for i=1:size(t1, 1)
        diffT = t1(i)-t2;
        gT(:, :, i) = 0.5 * D* kern.variance * (-h1(i,:).' + h2(:,i));
%                 + 2/(sigma*D*sqrt(pi))*exp(-diffT.*diffT/sigma2));
    end
end

if ~isfield(kern, 'isNormalised') || (kern.isNormalised == false)
    gT = gT * sigma * sqrt(pi);
end

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  error('simKerDiagGradX not implemented for gaussianInitial')
end
