function g = whiteblockKernGradient(kern, x, x2, covGrad)

% WHITEBLOCKKERNGRADIENT Gradient of WHITEBLOCK kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% white noise block kernel's parameters. As well as the kernel structure 
% and the input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
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
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO whiteblockKernParamInit, whiteblockKernDiagGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    covGrad = x2;
end

g = zeros(1, kern.nout);
if nargin < 4
    startOne = 1;
    endOne = 0;
    for i=1:kern.nout
        if iscell(x)
            endOne = endOne + size(x{i},1);
        else
            endOne = endOne + size(x,1);
        end
        g(1, i) = trace(covGrad(startOne:endOne, startOne:endOne));
        startOne = endOne + 1;
    end
end
