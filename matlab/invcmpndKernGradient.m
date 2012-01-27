function g = invcmpndKernGradient(kern, x, varargin)

% INVCMPNDKERNGRADIENT Gradient of INVERSE-PRESICION-CMPND kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% kernel's parameters theta_1, theta_2... . As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function L with respect to the
% relevant elements of the kernel matrix, i.e. dL/dK. 
% Then, dL/d theta = trace( dL/dK * dK/d theta)
% But dK/d theta_i = K ( inv(Ki) dKi/d theta_i inv(Ki) ) K, if theta_i is a
% parameter of Ki. We augment the partial derivatives matrix with the
% extra terms above, when we calculate dKi/ dtheta_i, to automatically get
% the correct expression from each function to xxxKernGradient, where xxx
% is the type of Ki.
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
% DESC computes the gradient of functions with respect to the
% kernel's parameters theta_1, theta_2... . As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function L with respect to the
% relevant elements of the kernel matrix, i.e. dL/dK. 
% Then, dL/d theta = trace( dL/dK * dK/d theta)
% But dK/d theta_i = K ( inv(Ki) dKi/d theta_i inv(Ki) ) K, if theta_i is a
% parameter of Ki. We augment the partial derivatives matrix with the
% extra terms above, when we calculate dKi/ dtheta_i, to automatically get
% the correct expression from each function to xxxKernGradient, where xxx
% is the type of Ki.
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
% SEEALSO invcmpndKernParamInit, kernGradient, invcmpndKernDiagGradient, kernGradX
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN

% !!! TODO: There are several precomputations and tricks to speed this up...!

g = zeros(1, size(kern.paramGroups, 1));
startVal = 1;
endVal = 0;

partial = varargin{end};
vararginOrig = varargin;

% Compute kernel. This might come from or be used as precomputation.
if nargin < 4
    K = kernCompute(kern, x);
else
    K = kernCompute(kern, x, varargin{1});
end
% !!! Maybe we can substitute some multiplications with .* !!
partial2 = K * partial * K; 


for i = 1:length(kern.comp)
    varargin = vararginOrig;
    endVal = endVal + kern.comp{i}.nParams;
     
    
    if ~isempty(kern.comp{i}.index)
        % only part of the data is involved in the kernel.
        if nargin < 4
            K_cur = kernCompute(kern.comp{i}, ...
                            x(:, kern.comp{i}.index));
            K_curInv = pdinv(K_cur);
            varargin{end} = K_curInv * partial2 * K_curInv;
            
            g(1, startVal:endVal)  = kernGradient(kern.comp{i}, ...
                         x(:, kern.comp{i}.index), ...
                         varargin{:});
        else
            K_cur = kernCompute(kern.comp{i}, ...
                         x(:, kern.comp{i}.index), ...
                         varargin{1}(:, kern.comp{i}.index));
            K_curInv = pdinv(K_cur);
            varargin{end} = K_curInv * partial2 * K_curInv;
            g(1, startVal:endVal) = kernGradient(kern.comp{i}, ...
                         x(:, kern.comp{i}.index), ...
                         varargin{1}(:, kern.comp{i}.index), ...
                         varargin{2:end});
        end
    else
        % all the data is involved with the kernel.
        if nargin < 4
            K_cur = kernCompute(kern.comp{i}, x);
        else
            K_cur = kernCompute(kern.comp{i}, x, varargin{1});
        end
        K_curInv = pdinv(K_cur);
        varargin{end} = K_curInv * partial2 * K_curInv;
        g(1, startVal:endVal)  = kernGradient(kern.comp{i}, x, varargin{:});
    end
    startVal = endVal + 1;
end
g = g*kern.paramGroups;

