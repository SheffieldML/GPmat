function [Ainv, UC] = pdinv(A, UC);

% PDINV Computes the inverse of a positive definite matrix

% IVM

if nargin < 2
  UC=[];
end
numData = size(A, 1);
try
  if isempty(UC)
    UC = chol(A);
  end
  invU = eye(numData)/UC;
  Ainv = invU*invU'; 
catch
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
  if nonPosDef
    warning(['Matrix is not positive definite in pdinv, adding jitter.'])
    
    lim = get(0, 'RecursionLimit');
    set(0, 'RecursionLimit', 10);
    try 
      [Ainv, UC] = pdinv(A+eye(size(A, 1))*1e-6);
    catch
      errMsg = lasterr;
      if findstr(errMsg, 'RecursionLimit')
	warning(['Jitter failed using SVD.'])
	[U, S, V] = svd(A);
	Ainv = V*diag(1./diag(S))*U';
	UC=[];
      else
	error(lasterr)
      end
    end
    set(0, 'RecursionLimit', lim);
    return
  else
    error(lasterr)
  end
end

