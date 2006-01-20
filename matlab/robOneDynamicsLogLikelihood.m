function ll = robOneDynamicsLogLikelihood(model)

% ROBONEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot one dynamics part.

% FGPLVM

thetaDiff = model.theta(2:end)-model.theta(1:end-1);
while any(thetaDiff>pi)
  ind = find(thetaDiff>pi);
  thetaDiff(ind) = thetaDiff(ind) - 2*pi;
end
while any(thetaDiff<-pi)
  ind = find(thetaDiff<-pi);
  thetaDiff(ind) = thetaDiff(ind) + 2*pi;
end

logLikeTheta = log(2*pi*model.sigma2)...
    +(thetaDiff.*thetaDiff)/model.sigma2;
likeTheta = exp(-0.5*logLikeTheta)*model.mixTheta;
likeTheta = likeTheta + (1-model.mixTheta)/(2*pi);
ll = sum(log(likeTheta));


logLikeR1 = model.a*log(model.b) - gammaln(model.a) ...
    +(model.a-1)*log(model.r) - model.b*model.r + log(model.mixR);
logLikeR2 = log(model.b) - model.b*model.r + log(1-model.mixR);

logLikeR1(find(logLikeR1<-316))=-316;
logLikeR2(find(logLikeR2<-316))=-316;

ind = find(logLikeR1>logLikeR2);
logLikeR(ind) = log(1+exp(logLikeR2(ind))./exp(logLikeR1(ind)))+logLikeR1(ind);

ind2 = find(logLikeR1<=logLikeR2);
logLikeR(ind2) = log(1+exp(logLikeR1(ind2))./exp(logLikeR2(ind2)))+logLikeR2(ind2);

ll = ll + sum(logLikeR);

