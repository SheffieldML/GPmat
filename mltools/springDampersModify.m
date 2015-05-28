function springDampers = springDampersModify(springDampers, springDamperVals, startPoints, ...
                                         staticPoints, widths, springConst, ...
                                         damperConst)

% SPRINGDAMPERSMODIFY Helper code for visualisation of springDamper data.
% FORMAT
% DESC is a helper function for visualising spring dampers using latent
% variable models.
% ARG springDampers : array of spring damper drawing objects.
% ARG springDamperValues : the values to set the springDamper data to.
% ARG startPoints : the fixed points of the spring damers.
% ARG widths : the widths of the spring dampers.
% ARG springConst : the spring constants of the spring dampers.
% ARG damperConst : the damper constants of the spring dampers.
% RETURN springDampers : contains array of spring damper objects.
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2006
%
% SEEALSO : springDamperVisualise, fgplvmResultsDynamic

% MLTOOLS


for i = 1:length(springDampers)
  vec = (staticPoints(i, :) - springDampers(i).start);
  vecNorm = vec./(sum(sqrt(vec.^2)));
  springDampers(i).end = springDamperVals(i)*vecNorm + staticPoints(i, :);
end
springDampers = springDamperDraw(springDampers);
