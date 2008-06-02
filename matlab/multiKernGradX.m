function gX_u = multiKernGradX(kern, x, x2, covGrad)

% MULTIKERNGRADX Gradient of MULTI kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the multiple output block
% kernel with respect to the input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x : locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in X.
%
% FORMAT
% DESC computes the gradident of the multiple output block
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in X1, numData2 is the number of data
% points in X2 and numInputs is the number of input
% dimensions in X.
%
% SEEALSO : multiKernParamInit, kernGradX, multiKernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008  

% KERN


if iscell(x)
  if nargin > 3 && ~iscell(x2)
    error('Time course information is not matched in Cell format!');
  end
  
  % Collate arguments.
  for i=1:kern.numBlocks
    dim1(i) = size(x{i}, 1);
    arg{i}{1} = x{i};
    if nargin > 3
      dim2(i) = size(x2{i}, 1);
      arg{i}{2} = x2{i};
    else
      dim2(i) = dim1(i);
      arg{i}{2} = arg{i}{1};
      covGrad = x2;
    end
  end
  
  gX = zeros(size( x{1},1),size( x{1},2),kern.numBlocks);

  cont = 1;
  for i = 1:kern.numBlocks

      startOne = sum(dim1(1:(i-1)))+1;
      endOne = sum(dim1(1:i));
      startThree = sum(dim2(1:(i-1))) + 1;
      endThree = sum(dim2(1:i));
      if i<= 1
          if nargin > 3
              gX(: ,: , cont) = multiKernGradientBlockX(kern, ...
                  arg{i}{:}, covGrad(startOne:endOne, ...
                  startThree:endThree), i, i);
          else
              gX(: ,: , cont) = multiKernGradientBlockX(kern, ...
                  arg{i}{1}, covGrad(startOne:endOne, ...
                  startThree:endThree), i, i);
          end
          cont = cont + 1;
      end
      for j = 1:i-1
          if j<=1
              if ~isempty(kern.block{i}.cross{j})
                  startTwo = sum(dim2(1:(j-1))) + 1;
                  endTwo =  sum(dim2(1:j));
                  gX(: , : , cont) = multiKernGradientBlockX(kern, arg{i}{1}, ...
                      arg{j}{2}, covGrad(startOne:endOne, ...
                      startTwo:endTwo), i, j);
                  cont = cont+1;
              end
          end
      end
  end
else
  % Something is wrong, because this is only allowed for sparse multigp
  
  % Collate arguments.
%   dim1 = size(x, 1);
%   arg{1} = x;
%   if nargin > 3
%     dim2 = size(x2, 1);
%     arg{2} = x2;
%   else
%     dim2 = dim1;
%     covGrad = x2;
%   end
%   
%   g = zeros(1, size(kern.paramGroups, 1));
%   startVal = 1;
%   endVal = 0;
%   for i = 1:kern.numBlocks
%     endVal = endVal + kern.comp{i}.nParams;
%     
%     startOne = (i-1)*dim1 + 1;
%     endOne = i*dim1;
%     if nargin > 3
%       g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
%                                                    arg{:}, ...
%                                                    covGrad(startOne:endOne, ...
%                                                       (i-1)*dim2 + 1:i*dim2), ...
%                                                  i, i);
%     else
%       g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
%                                                    arg{1}, ...
%                                                    covGrad(startOne:endOne, ...
%                                                       (i-1)*dim2 + 1:i*dim2), ...
%                                                  i, i);
%     end
%     startVal2 = 1;
%     endVal2 = 0;
%     for j = 1:i-1
%       endVal2 = endVal2 + kern.comp{j}.nParams;
%       if ~isempty(kern.block{i}.cross{j})
%         startTwo = (j-1)*dim2 + 1;
%         endTwo = j*dim2;
%         [g1, g2] = multiKernGradientBlock(kern, ...
%                                         arg{:}, ...
%                                         covGrad(startOne:endOne, ...
%                                                 startTwo:endTwo), ...
%                                         i, j);
% 
%         if nargin > 3
%           startThree = (j-1)*dim1 + 1;
%           endThree = j*dim1;
%           [g3 g4] = multiKernGradientBlock(kern, ...
%                                         arg{end:-1:1}, ...
%                                         covGrad(startThree:endThree, ...
%                                                 startTwo:endTwo)', j, i);
%           g(1, startVal:endVal) = g(1, startVal:endVal) + g1 + g4;
%           g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + g2 + g3;
%         else
%           g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
%           g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;           
%         end
%       end
%       startVal2 = endVal2 + 1;
%     end
%     startVal = endVal + 1;
%   end

end

gX_u = sum(gX, 3);
