function [K, P, Kbase,  dist, Factor, precColsInv, Pinv] = ...
    gaussianwhiteKernCompute(kern, x, x2)

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
    if kern.nIndFunct == 1,
        P = kern.precisionT/2;
        sqrtP = sqrt(kern.precisionT/2);        
        Px = x*sparseDiag(sqrtP);
        if nargin < 3
            n2 = dist2(Px, Px);
        else
            Px2 = x2*sparseDiag(sqrtP);
            n2 = dist2(Px, Px2);
        end
        Kbase = exp(-0.5*n2);
        K = kern.sigma2Noise*Kbase;
    else
        if nargin < 3
            x2 = x;
        end
        dist = zeros(size(x,1), size(x2,1), kern.inputDimension);
        k = zeros(size(x,1), size(x2,1), kern.inputDimension);
        kBase = zeros(size(x,1), size(x2,1), kern.inputDimension);
        if nargin < 3
            factor = zeros(size(x,1), size(x2,1), kern.inputDimension);
            P = zeros(size(x,1), size(x2,1), kern.inputDimension);
            Pinv = zeros(size(x,1), size(x2,1), kern.inputDimension);
            precColsInv = zeros(kern.nIndFunct, size(x,1), kern.inputDimension);
        end
        for i=1:kern.inputDimension,
            dist(:,:,i) = dist2(x(:,i), x2(:,i));
            if nargin < 3
                precCols = repmat(kern.precisionT(i,:)', 1, size(x,1));
                precRows = repmat(kern.precisionT(i,:) , size(x,1), 1);
                precColsInv(:,:,i) = 1./precCols;
                precRowsInv = 1./precRows;
                Pinv(:,:,i) = precColsInv(:,:,i) + precRowsInv;
                P(:,:,i) = 1./Pinv(:,:,i);
                precColInvNum =  repmat((kern.precisionT(i,:).^(-1/4))', 1, size(x,1));
                precRowsInvNum=  repmat(kern.precisionT(i,:).^(-1/4) , size(x,1), 1);
                factorDen = Pinv(:,:,i).^(1/2);
                factor(:,:,i) = 2^(1/2)*(precColInvNum.*precRowsInvNum)./factorDen;
                n2  = P(:,:,i).*dist(:,:,i);
                kBase(:,:,i) = exp(-0.5*n2);
                k(:,:,i) = factor(:,:,i).*kBase(:,:,i);
            else
                precCols = kern.precisionT(i,:)'/2;
                precColsInv = 1./precCols;
                Pinv =  precColsInv;
                P = repmat(1./Pinv', size(x,1), 1);
                % /~ MAURICIO : This is just for kernTest
                % P = repmat(1./Pinv, 1, size(x2,1));
                % ~/
                n2  = P.*dist(:,:,i);
                kBase(:,:,i) = exp(-0.5*n2);
                k(:,:,i) = kBase(:,:,i);
            end
        end
        K = k(:,:,1);
        Kbase = kBase(:,:,1);
        if nargin < 3
            Factor = factor(:,:,1);
        end
        for i=2:kern.inputDimension
            K = K.*k(:,:,i);
            Kbase = Kbase.*kBase(:,:,i);
            if nargin < 3
                Factor = Factor.*factor(:,:,i);
            end
        end
        K = kern.sigma2Noise*K;
    end
else
    if kern.nIndFunct == 1,
        if nargin<3
            x2 = x;
        end
        dist = dist2(x, x2);       
        P = kern.precisionT/2;
        Kbase = exp(-0.5.*P.*dist);
        K = kern.sigma2Noise*Kbase;        
    else
        if nargin < 3
            if kern.nIndFunct~=size(x,1)
                error(['The number of inducing functions must be equal the' ...
                    'number of inducing points']);
            end
            x2 = x;
        end
        dist = dist2(x, x2);
        if nargin < 3
            precCols = repmat(kern.precisionT', 1, size(x,1));
            precRows = repmat(kern.precisionT , size(x,1), 1);
            precColsInv = 1./precCols;
            precRowsInv = 1./precRows;
            Pinv = precColsInv + precRowsInv;
            P = 1./Pinv;
            precColInvNum =  repmat((kern.precisionT.^(-kern.inputDimension/4))', 1, size(x,1));
            precRowsInvNum=  repmat(kern.precisionT.^(-kern.inputDimension/4) , size(x,1), 1);
            factorDen = Pinv.^(kern.inputDimension/2);
            Factor = 2^(kern.inputDimension/2)*(precColInvNum.*precRowsInvNum)./factorDen;
            n2  = P.*dist;
            Kbase = exp(-0.5*n2);
            K = kern.sigma2Noise*Factor.*Kbase;
        else
            precCols = kern.precisionT'/2;
            precColsInv = 1./precCols;
            Pinv =  precColsInv;
            P = repmat(1./Pinv', size(x,1), 1);
            % /~ MAURICIO : This is just for kernTest
            % P = repmat(1./Pinv, 1, size(x2,1));
            % ~/
            n2  = P.*dist;
            Kbase = exp(-0.5*n2);
            K = kern.sigma2Noise*Kbase;
        end
    end
end
