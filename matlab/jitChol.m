function UC = jitChol(A)

% JITCHOL Do a Cholesky decomposition, if matrix isn't positive definite add jitter and do it again.

% NDLUTIL

maxTries = 10;
jitter = 0;
for i = 1:maxTries
  try
    % Try --- need to check A is positive definite
    if jitter == 0;
      jitter = abs(mean(diag(A)))*1e-6;
      UC = chol(A);
      break
    else
      warning(['Matrix is not positive definite in jitChol, adding ' num2str(jitter) ' jitter.'])
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
  else
    error(lasterr)
  end
end
