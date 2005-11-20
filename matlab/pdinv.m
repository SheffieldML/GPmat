function [Ainv, UC] = pdinv(A, UC);

% PDINV Computes the inverse of a positive definite matrix

% NDLUTIL

persistent JITTER;
maxTries = 10;
if nargin < 2
  UC=[];
end
jitter = 0;
% Try and invert using Cholesky decomposition.
if isempty(UC)
  % Cholesky wasn't passed, compute it.
  nonPosDef = 1;
  for i = 1:maxTries
    try
      % Try --- need to check A is positive definite
      if jitter == 0;
        jitter = mean(diag(A));
        UC = chol(A);
        break
      else
        warning(['Matrix is not positive definite in pdinv, adding ' num2str(jitter) ' jitter.'])
        UC = chol(A+jitter*eye(size(A, 1)));
        break
      end
    catch
      % Was the error due to not positive definite?
      nonPosDef = 0;
      verString = version;
      if str2double(verString(1:3)) > 6.1
        [void, errid] = lasterr;
        if strcmp(errid, 'MATLAB:posdef')
          nonPosDef = 1;
        end
      else
        errMsg = lasterr;
        if findstr(errMsg, 'positive definite')
          nonPosDef = 1;
        end
      end
    end
    if nonPosDef
      jitter = jitter*10;
    else
      error(lasterr)
    end
  end
end

invU = eye(size(A, 1))/UC;
Ainv = invU*invU'; 
