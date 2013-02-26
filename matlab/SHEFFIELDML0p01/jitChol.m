function [UC, jitter] = jitChol(A, maxTries)

% JITCHOL Do a Cholesky decomposition with jitter.
%
%	Description:
%
%	U = JITCHOL(A, MAXTRIES) attempts a Cholesky decomposition on the
%	given matrix, if matrix isn't positive definite the function gives a
%	warning, adds 'jitter' and tries again. At the first attempt the
%	amount of jitter added is 1e-6 times the mean of the diagonal.
%	Thereafter the amount of jitter is multiplied by 10 each time it is
%	added again. This is continued for a maximum of 10 times.
%	 Returns:
%	  U - the Cholesky decomposition for the matrix.
%	 Arguments:
%	  A - the matrix for which the Cholesky decomposition is required.
%	  MAXTRIES - the maximum number of times that jitter is added before
%	   giving up (default 10).
%
%	[U, JITTER] = JITCHOL(A, MAXTRIES) attempts a Cholesky decomposition
%	on the given matrix, if matrix isn't positive definite the function
%	adds 'jitter' and tries again. Thereafter the amount of jitter is
%	multiplied by 10 each time it is added again. This is continued for
%	a maximum of 10 times.  The amount of jitter added is returned.
%	 Returns:
%	  U - the Cholesky decomposition for the matrix.
%	  JITTER - the amount of jitter that was added to the matrix.
%	 Arguments:
%	  A - the matrix for which the Cholesky decomposition is required.
%	  MAXTRIES - the maximum number of times that jitter is added before
%	   giving up (default 10).
%	
%
%	See also
%	CHOL, PDINV, LOGDET


%	Copyright (c) 2005, 2006 Neil D. Lawrence


if nargin < 2
  maxTries = 10;
end
jitter = 0;
for i = 1:maxTries
  try
    % Try --- need to check A is positive definite
    if jitter == 0;
      jitter = abs(mean(diag(A)))*1e-6;
      UC = chol(A);
      break
    else
      if nargout < 2
        warning(['Matrix is not positive definite in jitChol, adding ' num2str(jitter) ' jitter.'])
      end
      UC = chol(real(A+jitter*eye(size(A, 1))));
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
    if i==maxTries
      error(['Matrix is non positive definite tried ' num2str(i) ...
             ' times adding jitter, but failed with jitter ' ...
             'of ' num2str(jitter) '. Increase max tries'])
    end
  else
    error(lasterr)
  end
end


