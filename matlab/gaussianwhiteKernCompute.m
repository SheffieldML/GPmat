function [K, Pinv]  = gaussianwhiteKernCompute(kern, x, x2)

% GAUSSIANWHITEKERNCOMPUTE Compute the covariance of the output samples 
% when the input is a white noise process and the smoothing kernel is a 
% Gaussian kernel.
% FORMAT
% DESC computes the kernel parameters for the Gaussian white kernel given
% inputs associated with rows and columns.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 : the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the Gaussian white kernel given a design matrix of inputs.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianwhiteKernParamInit, kernCompute, kernCreate, gaussianwhiteKernDiagCompute
% 
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if kern.isArd
    sqrtP = sqrt(kern.precisionT/2);
    Pinv = 2./kern.precisionT;
    Px = x*sparseDiag(sqrtP);
    if nargin < 3
        n2 = dist2(Px, Px);        
    else
        Px2 = x2*sparseDiag(sqrtP);
        n2 = dist2(Px, Px2);        
    end
    K = kern.sigma2Noise*exp(-0.5*n2);
else    
    if nargin < 3
        if kern.nIndFunct~=size(x,1)
            error(['The number of inducing functions must be equal the' ...
                'number of inducing points']);
        end
        x2 = x;
    end
    n2 = dist2(x, x2);
    if nargin < 3,
        precCols = repmat(kern.precisionT', 1, size(x,1));
        precRows = repmat(kern.precisionT , size(x,1), 1);
        precColsInv = 1./precCols;
        precRowsInv = 1./precRows;
        Pinv = precColsInv + precRowsInv;
        P = 1./Pinv;
        precColInvNum =  repmat((kern.precisionT.^(-kern.inputDimension/4))', 1, size(x,1));
        precRowsInvNum=  repmat(kern.precisionT.^(-kern.inputDimension/4) , size(x,1), 1);
        factorDen = Pinv.^(kern.inputDimension/2); 
        factor = 2^(kern.inputDimension/2)*(precColInvNum.*precRowsInvNum)./factorDen;
        K = kern.sigma2Noise.*factor.*exp(-0.5.*P.*n2);
    else
        precCols = kern.precisionT'/2;
        precColsInv = 1./precCols;
        Pinv =  precColsInv;
        P = repmat(1./Pinv', size(x,1), 1);    
        %P = repmat(1./Pinv, 1, size(x2,1));
        K = kern.sigma2Noise*exp(-0.5.*P.*n2);        
    end
    
end
