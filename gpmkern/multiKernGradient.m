function g = multiKernGradient(kern, x, x2, covGrad)

% MULTIKERNGRADIENT Gradient of MULTI kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% multiple output block
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
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
% SEEALSO : multiKernParamInit, kernGradient, multiKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007
  
% KERN

if iscell(x)
  if nargin > 3 && ~iscell(x2)
    error('Time course information is not matched in Cell format!');
  end
  
  % Collate arguments.
  for i=1:kern.numBlocks
    % The vector dim1 contains the row size of each block,
    % and arg{i}{1} contains the row inputs of the ith block
    dim1(i) = size(x{i}, 1);
    arg{i}{1} = x{i};
    if nargin > 3
      % The vector dim2 contains the column size of each block,
      % and arg{i}{1} contains the column inputs of the ith block
      dim2(i) = size(x2{i}, 1);
      arg{i}{2} = x2{i};
    else
      % If separate column inputs are not given they are the same
      % as the row inputs
      dim2(i) = dim1(i);
      arg{i}{2} = arg{i}{1};

      % If the column inputs were not given, then covGrad was given in the x2 parameter
      covGrad = x2;
    end
  end

  %x
  %x2
  %dim1
  %dim2
  
  g = zeros(1, size(kern.paramGroups, 1));
  startVal = 1;
  endVal = 0;
  
  for i = 1:kern.numBlocks
    % Compute contribution of the ith diagonal block to the gradient
    endVal = endVal + kern.comp{i}.nParams;
    
    startOne = sum(dim1(1:(i-1)))+1;
    endOne = sum(dim1(1:i));
    startThree = sum(dim2(1:(i-1))) + 1;
    endThree = sum(dim2(1:i));

    % Compute the contribution of the block if it has nonzero rows and columns
    if ((endOne >= startOne) && (endThree >= startThree)),
      if nargin > 3  % If x2 has been provided use those indices
        g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
                              arg{i}{:}, covGrad(startOne:endOne, ...
                              startThree:endThree), i, i);
      else
        g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
                              arg{i}{1}, covGrad(startOne:endOne, ...
                              startThree:endThree), i, i);
      end
    else
      % If the block has zero rows or columns then its gradient is all zero
      g(1, startVal:endVal) = 0;
    end;

    % Compute gradient contributions of the off-diagonal blocks in the
    % same row/column as the ith diagonal block
    startVal2 = 1;
    endVal2 = 0;
    for j = 1:i-1
      endVal2 = endVal2 + kern.comp{j}.nParams;
      if ~isempty(kern.block{i}.cross{j})
        startTwo = sum(dim2(1:(j-1))) + 1;
        endTwo =  sum(dim2(1:j));

        % Compute the contribution of the block if it has nonzero rows and columns
        if ((endOne >= startOne) && (endTwo >= startTwo)),
          [g1, g2] = multiKernGradientBlock(kern, arg{i}{1}, ...
                              arg{j}{2}, covGrad(startOne:endOne, ...
                              startTwo:endTwo), i, j);
        
          g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
          g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;
        end;
      end
      startVal2 = endVal2 + 1;
    end
    startVal = endVal + 1;
  end 
  
else
  
  % Collate arguments.
  dim1 = size(x, 1);
  arg{1} = x;
  if nargin > 3
    dim2 = size(x2, 1);
    arg{2} = x2;
  else
    dim2 = dim1;
    covGrad = x2;
  end
  
  g = zeros(1, size(kern.paramGroups, 1));
  startVal = 1;
  endVal = 0;
  for i = 1:kern.numBlocks
    endVal = endVal + kern.comp{i}.nParams;
    
    startOne = (i-1)*dim1 + 1;
    endOne = i*dim1;
    if nargin > 3
      g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
                                                   arg{:}, ...
                                                   covGrad(startOne:endOne, ...
                                                      (i-1)*dim2 + 1:i*dim2), ...
                                                 i, i);
    else
      g(1, startVal:endVal) = multiKernGradientBlock(kern, ...
                                                   arg{1}, ...
                                                   covGrad(startOne:endOne, ...
                                                      (i-1)*dim2 + 1:i*dim2), ...
                                                 i, i);
    end
    startVal2 = 1;
    endVal2 = 0;
    for j = 1:i-1
      endVal2 = endVal2 + kern.comp{j}.nParams;
      if ~isempty(kern.block{i}.cross{j})
        startTwo = (j-1)*dim2 + 1;
        endTwo = j*dim2;
        [g1, g2] = multiKernGradientBlock(kern, ...
                                        arg{:}, ...
                                        covGrad(startOne:endOne, ...
                                                startTwo:endTwo), ...
                                        i, j);

        if nargin > 3
          startThree = (j-1)*dim1 + 1;
          endThree = j*dim1;
          [g3 g4] = multiKernGradientBlock(kern, ...
                                        arg{end:-1:1}, ...
                                        covGrad(startThree:endThree, ...
                                                startTwo:endTwo)', j, i);
          g(1, startVal:endVal) = g(1, startVal:endVal) + g1 + g4;
          g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + g2 + g3;
        else
          g(1, startVal:endVal) = g(1, startVal:endVal) + 2*g1;
          g(1, startVal2:endVal2) = g(1, startVal2:endVal2) + 2*g2;           
        end
      end
      startVal2 = endVal2 + 1;
    end
    startVal = endVal + 1;
  end

end

g = g*kern.paramGroups;
