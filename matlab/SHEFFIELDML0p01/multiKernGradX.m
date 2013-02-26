function gX = multiKernGradX(kern, x, x2)

% MULTIKERNGRADX Gradient of MULTI kernel with respect to a point x.
%
%	Description:
%
%	GX = MULTIKERNGRADX(KERN, X) computes the gradient of the multiple
%	output block kernel with respect to the input positions.
%	 Returns:
%	  GX - the returned gradients. The gradients are returned in a
%	   matrix which is numData x numInputs x numData. Where numData is
%	   the number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = MULTIKERNGRADX(KERN, X, X2) computes the gradident of the
%	multiple output block kernel with respect to the input positions
%	where both the row positions and column positions are provided
%	separately.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData2 x numInputs x numData1. Where numData1 is the
%	   number of data points in X1, numData2 is the number of data points
%	   in X2 and numInputs is the number of input dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, KERNGRADX, MULTIKERNDIAGGRADX


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


if iscell(x)
    if nargin > 3 && ~iscell(x2)
        error('Time course information is not matched in Cell format!');
    end
    % Collate arguments.
    for i=1:kern.numBlocks
        arg{i}{1} = x{i};
        if nargin > 2
            arg{i}{2} = x2{i};
        else
            arg{i}{2} = arg{i}{1};
        end
    end
    for i = 1:kern.numBlocks
        if nargin > 2
            gX = multiKernGradientBlockX(kern, arg{i}{:}, i, i);
        else
            gX = multiKernGradientBlockX(kern, arg{i}{1}, i, i);
        end
        for j = 1:i-1
            if ~isempty(kern.block{i}.cross{j})
                gX = multiKernGradientBlockX(kern, arg{i}{1}, arg{j}{2}, i, j);
            end
        end
    end
end

