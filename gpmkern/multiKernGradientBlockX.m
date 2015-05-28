function gX = multiKernGradientBlockX(kern, X, X2, covGrad, i, j)

% MULTIKERNGRADIENTBLOCKX
% FORMAT
% DESC computes the gradient with respect to the inducing points for a block of a
% multi-output kerneal given two matrices of input.
% ARG kern : the structure containing the kernel.
% ARG X1 : first set of kernel inputs.
% ARG X2 : second set of kernel inputs.
% ARG covGrad : Gradient of the objective function with respect to
% the relevant portion of the kernel matrix.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% RETURN g1 : the gradient of the kernel parameters from the first
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN g2 : the gradient of the kernel parameters from the second
% kernel in the order provided by the relevant kernExtractParam commands.
%
% FORMAT
% DESC compute a block of a multi-output kernel given a single matrix
% of input.
% ARG kern : the structure containing the kernel.
% ARG X : first set of kernel inputs.
% ARG covGrad : Gradient of the objective function with respect to
% the relevant portion of the kernel matrix.
% ARG i : the row of the block of the kernel to be computed.
% ARG j : the column of the block of the kernel to be computed.
% RETURN g1 : the gradient of the kernel parameters from the first
% kernel in the order provided by the relevant kernExtractParam commands.
% RETURN g2 : the gradient of the kernel parameters from the second
% kernel in the order provided by the relevant kernExtractParam commands.
%
% SEEALSO : multiKernCreate, multiKernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

% KERN

if nargin < 6
  j = i;
  i = covGrad;
  covGrad = X2;
  X2 = [];
end

if i == j
    fhandle = [kern.comp{i}.type 'KernGradX'];
    arg{1} = kern.comp{i};
else
    if j<i
        fhandle = [kern.block{i}.cross{j} 'KernGradX'];
        transpose = kern.block{i}.transpose(j);
    else
        fhandle = [kern.block{j}.cross{i} 'KernGradX'];
        transpose = ~kern.block{j}.transpose(i);
    end
    if transpose
        arg{1} = kern.comp{j};
        arg{2} = kern.comp{i};
    else
        arg{1} = kern.comp{i};
        arg{2} = kern.comp{j};
    end
end
fhandle = str2func(fhandle);
arg{end+1} = X;
if ~isempty(X2);
    arg{end+1} = X2;
end

gX = fhandle(arg{:}, covGrad);

%/~ the old version of the code
% if nargin < 6
%   j = i;
%   i = covGrad;
%   covGrad = X2;
%   X2 = [];
% end
% 
% outArg = 2;
% if i == j
%     fhandle = [kern.comp{i}.type 'KernGradX'];
%     transpose = 0;
%     arg{1} = kern.comp{i};
%     factors = kernFactors(kern.comp{i}, 'gradfact');
%     outArg = 1;
% else
%     if j<i
%         fhandle = [kern.block{i}.cross{j} 'KernGradX'];
%         transpose = kern.block{i}.transpose(j);
%     else
%         fhandle = [kern.block{j}.cross{i} 'KernGradX'];
%         transpose = ~kern.block{j}.transpose(i);
%     end
%     if transpose
%         arg{1} = kern.comp{j};
%         factors{1} = kernFactors(kern.comp{j}, 'gradfact');
%         arg{2} = kern.comp{i};
%         factors{2} = kernFactors(kern.comp{i}, 'gradfact');
%     else
%         arg{1} = kern.comp{i};
%         factors{1} = kernFactors(kern.comp{i}, 'gradfact');
%         arg{2} = kern.comp{j};
%         factors{2} = kernFactors(kern.comp{j}, 'gradfact');
%     end
% end
% fhandle = str2func(fhandle);
% arg{end+1} = X;
% if ~isempty(X2);
%     arg{end+1} = X2;
% end
% switch outArg
%     case 1
%         gX = fhandle(arg{:}, covGrad);
%     case 2
%         gX = fhandle(arg{:}, covGrad);
%         if transpose
%             g = g2;
%             g2 = g1;
%             g1 = g;
%         end
%     otherwise
%         error('Invalid number of out arguments.')
% end
%~/
