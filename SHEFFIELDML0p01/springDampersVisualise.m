function springDampers = springDampersVisualise(springDamperVals, startPoints, ...
                                         staticPoints, widths, springConst, ...
                                         damperConst)

% SPRINGDAMPERSVISUALISE Helper code for showing an spring dampers during 2-D visualisation.
%
%	Description:
%
%	SPRINGDAMPERS = SPRINGDAMPERSVISUALISE(SPRINGDAMPERVALUES,
%	STARTPOINTS, WIDTHS, SPRINGCONST, DAMPERCONST) is a helper function
%	for plotting spring damper locations using latent variable models.
%	 Returns:
%	  SPRINGDAMPERS - contains array of spring damper objects.
%	 Arguments:
%	  SPRINGDAMPERVALUES - the values to set the springDamper data to.
%	  STARTPOINTS - the fixed points of the spring damers.
%	  WIDTHS - the widths of the spring dampers.
%	  SPRINGCONST - the spring constants of the spring dampers.
%	  DAMPERCONST - the damper constants of the spring dampers.
%	
%
%	See also
%	SPRINGDAMPERSMODIFY, FGPLVMRESULTSDYNAMIC


%	Copyright (c) 2007 Neil D. Lawrence

  
for i = 1:length(springDamperVals)
  springDampers(i) = springDamperCreate(startPoints(i, :), ...
                                        staticPoints(i, :), ...
                                        widths(i), springConst(i), ...
                                        damperConst(i));
  springDampers(i).selected = 0;
end
springDampers = springDampersModify(springDampers, springDamperVals, ...
                                       startPoints, staticPoints, widths, ...
                                       springConst, damperConst);
for i = 1:length(springDamperVals)
  set(springDampers(i).handle, 'linewidth', 2)
  set(springDampers(i).spring.handle, 'linewidth', 2)
  set(springDampers(i).damper.handle, 'linewidth', 2)
end