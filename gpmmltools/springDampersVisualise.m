function springDampers = springDampersVisualise(springDamperVals, startPoints, ...
                                         staticPoints, widths, springConst, ...
                                         damperConst)

% SPRINGDAMPERSVISUALISE Helper code for showing an spring dampers during 2-D visualisation.
% FORMAT
% DESC is a helper function for plotting spring damper locations using latent
% variable models.
% ARG springDamperValues : the values to set the springDamper data to.
% ARG startPoints : the fixed points of the spring damers.
% ARG widths : the widths of the spring dampers.
% ARG springConst : the spring constants of the spring dampers.
% ARG damperConst : the damper constants of the spring dampers.
% RETURN springDampers : contains array of spring damper objects.
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% SEEALSO : springDampersModify, fgplvmResultsDynamic

% MLTOOLS
  
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
