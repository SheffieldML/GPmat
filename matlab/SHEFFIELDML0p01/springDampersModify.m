function springDampers = springDampersModify(springDampers, springDamperVals, startPoints, ...
                                         staticPoints, widths, springConst, ...
                                         damperConst)

% SPRINGDAMPERSMODIFY Helper code for visualisation of springDamper data.
%
%	Description:
%
%	SPRINGDAMPERS = SPRINGDAMPERSMODIFY(SPRINGDAMPERS,
%	SPRINGDAMPERVALUES, STARTPOINTS, WIDTHS, SPRINGCONST, DAMPERCONST)
%	is a helper function for visualising spring dampers using latent
%	variable models.
%	 Returns:
%	  SPRINGDAMPERS - contains array of spring damper objects.
%	 Arguments:
%	  SPRINGDAMPERS - array of spring damper drawing objects.
%	  SPRINGDAMPERVALUES - the values to set the springDamper data to.
%	  STARTPOINTS - the fixed points of the spring damers.
%	  WIDTHS - the widths of the spring dampers.
%	  SPRINGCONST - the spring constants of the spring dampers.
%	  DAMPERCONST - the damper constants of the spring dampers.
%	
%
%	See also
%	SPRINGDAMPERVISUALISE, FGPLVMRESULTSDYNAMIC


%	Copyright (c) 2003, 2004, 2006 Neil D. Lawrence



for i = 1:length(springDampers)
  vec = (staticPoints(i, :) - springDampers(i).start);
  vecNorm = vec./(sum(sqrt(vec.^2)));
  springDampers(i).end = springDamperVals(i)*vecNorm + staticPoints(i, :);
end
springDampers = springDamperDraw(springDampers);
