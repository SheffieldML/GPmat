function [ld, U] = logdet(A, U)

% LOGDET The log of the determinant when argument is positive definite.

% IVM

if nargin < 2
  U = [];
end
try
  if isempty(U);
    U = chol(A);
  end
  ld = 2*sum(log(diag(U)));
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
    warning(['Matrix is not positive definite in pdinv, adding jitter'])
    lim = get(0, 'RecursionLimit');
    set(0, 'RecursionLimit', 10);
    try 
      [ld, U] = logdet(A+eye(size(A, 1))*1e-6);
    catch
      errMsg = lasterr;
      if findstr(errMsg, 'RecursionLimit')
	warning(['Jitter failed using log(det()).'])
	ld=log(det(A));
	U = [];
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

