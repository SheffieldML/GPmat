function K = multiKernCompute(kern, varargin)
% MULTIKERNCOMPUTE Compute the MULTI kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the multiple output block
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the multiple output block
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : multiKernParamInit, kernCompute, kernCreate, multiKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007
  
% KERN

if iscell(varargin{1})  % If x is a cell array of input points
  if length(varargin{1}) ~= kern.numBlocks
    error('Time information is not matched among blocks!');
  end
  for i = 1:kern.numBlocks
    dim1(i) = size(varargin{1}{i}, 1);    
    if length(varargin)>1  % If x2 is provided it must also be a cell array
      if length(varargin{1}) ~= length(varargin{2})
        error('Time information is not matched within the block!');
      end      
      dim2(i) = size(varargin{2}{i}, 1);
    else
      dim2(i) = dim1(i); % If x2 is provided is not provided just use x1
    end    
    
  end
    
  K = zeros(sum(dim1), sum(dim2));
  
  for i = 1:kern.numBlocks
    % Compute the i:th diagonal block
    startOne = sum(dim1(1:(i-1))) + 1;
    endOne = sum(dim1(1:i));
    startThree = sum(dim2(1:(i-1))) + 1;
    endThree = sum(dim2(1:i));

    % Compute the block if it has a nonzero amount of rows and columns
    if ((endOne >= startOne) && (endThree >= startThree)),
      if length(varargin)<2
        K(startOne:endOne, startThree:endThree) = multiKernComputeBlock(kern, ...
                                                      varargin{1}{i}, i, i);     
      else
        K(startOne:endOne, startThree:endThree) = multiKernComputeBlock(kern, ...
                                       varargin{1}{i}, varargin{2}{i}, i, i);      
      end
    end;

    % Compute the off-diagonal blocks in the same row/column as the i:th diagonal block
    for j = 1:i-1
      if ~isempty(kern.block{i}.cross{j})
        startTwo = sum(dim2(1:(j-1))) + 1;
        endTwo =  sum(dim2(1:j));
        
        if length(varargin)<2
          % If a separate set of column-inputs has not been given, 
          % row and column inputs are the same

	  %size(K)
	  %startOne
	  %endOne
	  %startTwo
	  %endTwo
	  %size(varargin{1})
	  %i
	  %j
	  %K(startOne:endOne, startTwo:endTwo)
	  %varargin{1}{i}
	  %varargin{1}{j}
          %kern
          %tempK = multiKernComputeBlock(kern, ...
          %                      varargin{1}{i}, varargin{1}{j}, i, j)

	  % Compute the block if it has a nonzero amount of rows and columns
	  if ((endOne >= startOne) && (endTwo >= startTwo)),
            K(startOne:endOne, startTwo:endTwo) = multiKernComputeBlock(kern, ...
                                  varargin{1}{i}, varargin{1}{j}, i, j);
  	    % Since the row inputs are the same as the column inputs, the 
            % corresponding block on the other side of the diagonal is just 
            % the transpose of the one just computed.
            K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                                                  startTwo:endTwo)';
          end;
        else
          % If a separate set of column-inputs has been given, use that

	  % Compute the block if it has a nonzero amount of rows and columns
	  if ((endOne >= startOne) && (endTwo >= startTwo)),
            K(startOne:endOne, startTwo:endTwo) = multiKernComputeBlock(kern, ...
                                  varargin{1}{i}, varargin{2}{j}, i, j);
          end;

	  % Since the row inputs are not the same as the column inputs, the 
          % corresponding block on the other side of the diagonal is not just
          % a transpose, so we compute it separately here.
          startFour = sum(dim1(1:(j-1))) + 1;
          endFour =  sum(dim1(1:j));
	  % Compute the block if it has a nonzero amount of rows and columns
	  if ((endFour >= startFour) && (endThree >= startThree)),
            K(startFour:endFour, startThree:endThree) = ...
                multiKernComputeBlock(kern, varargin{2}{i}, varargin{1}{j}, j, i)';
          end;
        end
      end
    end
  end
 
else
  dim1 = size(varargin{1}, 1);
  
  if length(varargin)>1
    dim2 = size(varargin{2}, 1);
  else
    dim2 = dim1;
  end
  
  K = zeros(kern.numBlocks*dim1, kern.numBlocks*dim2);
  
  for i = 1:kern.numBlocks
    startOne = (i-1)*dim1 + 1;
    endOne = i*dim1;
    startThree = (i-1)*dim2 + 1;
    endThree = i*dim2;
    K(startOne:endOne, startThree:endThree) = multiKernComputeBlock(kern, ...
                                                      varargin{:}, i, i);
    for j = 1:i-1
      if ~isempty(kern.block{i}.cross{j})
        startTwo = (j-1)*dim2 + 1;
        endTwo = j*dim2;
        K(startOne:endOne, startTwo:endTwo) = multiKernComputeBlock(kern, ...
                                                          varargin{:}, i, j);
        if length(varargin)<2
          K(startTwo:endTwo, startOne:endOne) = K(startOne:endOne, ...
                                                startTwo:endTwo)';
        else
          startFour = (j-1)*dim1 + 1;
          endFour = j*dim1;
          K(startFour:endFour, startThree:endThree) = ...
              multiKernComputeBlock(kern, varargin{end:-1:1}, j, i)';
        end
      end
    end
  end
  
end
