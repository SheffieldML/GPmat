function str = optimiDefaultConstraint(constraint)

% OPTIMIDEFAULTCONSTRAINT Returns function for parameter constraint.
% FORMAT
% DESC returns the current default function for constraining a
% parameter. Formerly (up to version OPTIMI 0.163) this was
% 'negLogLogit' for positive constraints, as this keeps things roughly linear in the
% positive half space, however, it is more standard to use 'exp'
% (i.e. optimise in the log space). This function allows you to
% change the option globally. The defaults are given below.
% ARG constraint : the type of constraint you want to place on the
% parameter, options include 'positive' (gives an 'exp' constraint)
% and 'zeroone' (gives a 'sigmoid' constraint).
% RETURN str : the type of function used to apply the constraint
% from the 'optimi' toolbox.
%
% SEEALSO : expTransform, sigmoidTransform, linearTransform, negLogLogitTransform
% COPYRIGHT : Neil D. Lawrence, 2006

% OPTIMI

switch constraint
 case 'positive'
  str = 'exp';
 case 'zeroone'
  str = 'sigmoid';
 case 'bounded'
  str = 'sigmoidab';
end
