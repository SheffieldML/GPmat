function [x, options, flog, pointlog, scalelog] = scg2(f, x, options, gradf, varargin)
% SCG2 Scaled conjugate gradient optimization like netlab's scg, with a slight modification for speeding it up.
% FORMAT
% DESC This is a slight modification of the netlab function for scg (from Ian T Nabney) so
% that it can (optionally) receive an ObjectiveGradient function (which returns both,
% the scalar value for the objective function and the vector of partial derivatives)
% for the first argument 'f', instead of an objective function. This
% means that, when possible, the precomputations will be made only once and used for
% both, the objective and the gradient function. However, the gradient function must also
% be provided, because in a certain step of the algorithm scg is interested in just 
% finding the gradient(whereas every calculation of the objective is always
% accompanyied by a calculation of the gradient). 
% This function can work as the traditional scg if 'f' is a function with only one
% output argument, otherwise, (if 'f' is a function with more than one outputs) 
% the behaviour described above will take effect. This is done, therefore, only
% in the places of the original scg where f and gradf are evaluated with the same
% set of parameters x (before they change in the next optimisation cycle).
% Details:
%
%	[X, OPTIONS] = SCG2(F, X, OPTIONS, GRADF) uses a scaled conjugate
%	gradients algorithm to find a local minimum of the a function
%	whose gradient is given by GRADF(X), where X is a row vector. The
%   function for which a minimum is found can either be F = F(X) (when the
%   provided F returns a scalar value) or the provided F can return both,
%   the scalar value for the objective function and the partial derivatives (as
%   GRADF). In the later case, some speed-up can be achieved.
%	The point at which the function has a local minimum is
%	returned as X.  The function value at that point is returned in
%	OPTIONS(8).
%
%	[X, OPTIONS, FLOG, POINTLOG, SCALELOG] = SCG2(F, X, OPTIONS, GRADF)
%	also returns (optionally) a log of the function values after each
%	cycle in FLOG, a log of the points visited in POINTLOG, and a log of
%	the scale values in the algorithm in SCALELOG.
%
%	SCG2(F, X, OPTIONS, GRADF, P1, P2, ...) allows additional arguments to
%	be passed to F() and GRADF().     The optional parameters have the
%	following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG, and the points visited in the
%	return argument POINTSLOG.  If OPTIONS(1) is set to 0, then only
%	warning messages are displayed.  If OPTIONS(1) is -1, then nothing is
%	displayed.
%
%	OPTIONS(2) is a measure of the absolute precision required for the
%	value of X at the solution.  If the absolute difference between the
%	values of X between two successive steps is less than OPTIONS(2),
%	then this condition is satisfied.
%
%	OPTIONS(3) is a measure of the precision required of the objective
%	function at the solution.  If the absolute difference between the
%	objective function values between two successive steps is less than
%	OPTIONS(3), then this condition is satisfied. Both this and the
%	previous condition must be satisfied for termination.
%
%	OPTIONS(9) is set to 1 to check the user defined gradient function.
%
%	OPTIONS(10) returns the total number of function evaluations
%	(including those in any line searches).
%
%	OPTIONS(11) returns the total number of gradient evaluations.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	See also
%	SCG, CONJGRAD, QUASINEW
%
%	Copyright (c) Ian T Nabney (1996-2001)
%   Modifications: Andreas C. Damianou 2011

% OPTIMI



%  Set up the options.
if length(options) < 18
  error('Options vector too short')
end

if(options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
gradcheck = options(9);

% Set up strings for evaluating function and gradient
f = fcnchk(f, length(varargin)); % f = @vargplvmObjective
gradf = fcnchk(gradf, length(varargin)); % gradf = @vargplvmGradient

nparams = length(x);

%  Check gradients
if (gradcheck)
  feval('gradchek', x, f, gradf, varargin{:});
end

sigma0 = 1.0e-4;
if  nargout(f) > 1
   [fold, gradnew] = feval(f, x, varargin{:}); 
else
    fold = feval(f, x, varargin{:});	% Initial function value.
    gradnew = feval(gradf, x, varargin{:});	% Initial gradient.
end
fnow = fold;
options(10) = options(10) + 1;		% Increment function evaluation counter.
%gradnew = feval(gradf, x, varargin{:});	% Initial gradient. 
gradold = gradnew;
options(11) = options(11) + 1;		% Increment gradient evaluation counter.
d = -gradnew;				% Initial search direction.
success = 1;				% Force calculation of directional derivs.
nsuccess = 0;				% nsuccess counts number of successes.
beta = 1.0;				% Initial scale parameter.
betamin = 1.0e-15; 			% Lower bound on scale.
betamax = 1.0e100;			% Upper bound on scale.
j = 1;					% j counts number of iterations.
if nargout >= 3
  flog(j, :) = fold;
  if nargout == 4
    pointlog(j, :) = x;
  end
end

% Main optimization loop.
while (j <= niters)

  % Calculate first and second directional derivatives.
  if (success == 1)
    mu = d*gradnew';
    if (mu >= 0)
      d = - gradnew;
      mu = d*gradnew';
    end
    kappa = d*d';
    if kappa < eps
      options(8) = fnow;
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    
    gplus = feval(gradf, xplus, varargin{:}); %Evaluate gradf
    options(11) = options(11) + 1; 
    theta = (d*(gplus' - gradnew'))/sigma;
  end

  % Increase effective curvature and evaluate step size alpha.
  delta = theta + beta*kappa;
  if (delta <= 0) 
    delta = beta*kappa;
    beta = beta - theta/kappa;
  end
  alpha = - mu/delta;
  
  % Calculate the comparison ratio.
  xnew = x + alpha*d;
  if nargout(f) > 1
    % (Modification): This causes the objectiveGradient to be called
    [fnew, gradf2] = feval(f, xnew, varargin{:});
  else
    fnew = feval(f, xnew, varargin{:}); %Evaluate f (objective).
  end
  options(10) = options(10) + 1;
  Delta = 2*(fnew - fold)/(alpha*mu);
  if (Delta  >= 0)
    success = 1;
    nsuccess = nsuccess + 1;
    x = xnew;
    fnow = fnew;
  else
    success = 0;
    fnow = fold;
  end

  if nargout >= 3
    % Store relevant variables
    flog(j) = fnow;		% Current function value
    if nargout >= 4
      pointlog(j,:) = x;	% Current position
      if nargout >= 5
	scalelog(j) = beta;	% Current scale parameter
      end
    end
  end    
  if display > 0
    fprintf(1, 'Cycle %4d  Error %11.6f  Scale %e\n', j, fnow, beta);
  end
  
  
  % (Modification): If success==1 (in the vast majority of iterations, normally),
  % then from above, x=xnew, so in this case the
  % calculation for F and G need to be done for the same set of params and
  % in that case we can calculate F and GRADF with the same
  % precomputations, so an objectiveGradient instead of an objective can
  % speed-up the algorithm at this point.
  if (success == 1)
    % Test for termination

    if (max(abs(alpha*d)) < options(2) & max(abs(fnew-fold)) < options(3))
      options(8) = fnew;
      return;
    else
      % Update variables for new position
      fold = fnew;
      gradold = gradnew;
      if nargout(f) > 1
          gradnew = gradf2; %objectiveGradient was used, so the gradient is already calculated
      else
          gradnew = feval(gradf, x, varargin{:}); % Calculate gradf
      end
      options(11) = options(11) + 1;
      % If the gradient is zero then we are done.
      if (gradnew*gradnew' == 0)
           options(8) = fnew;
           return;
      end
    end
  end

  % Adjust beta according to comparison ratio.
  if (Delta < 0.25)
    beta = min(4.0*beta, betamax);
  end
  if (Delta > 0.75)
    beta = max(0.5*beta, betamin);
  end

  % Update search direction using Polak-Ribiere formula, or re-start 
  % in direction of negative gradient after nparams steps.
  if (nsuccess == nparams)
    d = -gradnew;
    nsuccess = 0;
  else
    if (success == 1)
      gamma = (gradold - gradnew)*gradnew'/(mu);
      d = gamma*d - gradnew;
    end
  end
  j = j + 1;
end

% If we get here, then we haven't terminated in the given number of 
% iterations.

options(8) = fold;
if (options(1) >= 0)
  disp(maxitmess);
end

