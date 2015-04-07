function g = gpsimMarginalGradients(param, model)

model = gpsimMapExpandParam(model, param);

g = gpsimMapLogLikeGradients(model);