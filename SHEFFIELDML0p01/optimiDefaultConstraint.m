function str = optimiDefaultConstraint(constraint)

% OPTIMIDEFAULTCONSTRAINT Returns function for parameter constraint.
%
%	Description:
%
%	STR = OPTIMIDEFAULTCONSTRAINT(CONSTRAINT) returns the current
%	default function for constraining a parameter. Formerly (up to
%	version OPTIMI 0.163) this was 'negLogLogit' for positive
%	constraints, as this keeps things roughly linear in the positive
%	half space, however, it is more standard to use 'exp' (i.e. optimise
%	in the log space). This function allows you to change the option
%	globally. The defaults are given below.
%	 Returns:
%	  STR - the type of function used to apply the constraint from the
%	   'optimi' toolbox.
%	 Arguments:
%	  CONSTRAINT - the type of constraint you want to place on the
%	   parameter, options include 'positive' (gives an 'exp' constraint)
%	   and 'zeroone' (gives a 'sigmoid' constraint).
%
%	See also
%	EXPTRANSFORM, SIGMOIDTRANSFORM, LINEARTRANSFORM, NEGLOGLOGITTRANSFORM


%	Copyright (c) 2006 Neil D. Lawrence


switch constraint
 case 'positive'
  str = 'exp';
 case 'zeroone'
  str = 'sigmoid';
 case 'bounded'
  str = 'sigmoidab';
end
