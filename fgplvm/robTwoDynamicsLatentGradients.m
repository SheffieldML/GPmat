function gX = robTwoDynamicsLatentGradients(model);

% ROBTWODYNAMICSLATENTGRADIENTS Gradients of the X vector given the dynamics model.

% FGPLVM

gX = zeros(size(model.X));

thetaDiff = model.theta(2:end)-model.theta(1:end-1);
while any(thetaDiff>pi)
  ind = find(thetaDiff>pi);
  thetaDiff(ind) = thetaDiff(ind) - 2*pi;
end
while any(thetaDiff<-pi)
  ind = find(thetaDiff<-pi);
  thetaDiff(ind) = thetaDiff(ind) + 2*pi;
end

logLikeThetaOne = log(2*pi*model.sigma2One)...
    +(thetaDiff.*thetaDiff)/model.sigma2One;
gLikeThetaOne = exp(-0.5*logLikeThetaOne)*model.mixThetaOne;
logLikeThetaTwo = log(2*pi*model.sigma2Two)...
    +(thetaDiff.*thetaDiff)/model.sigma2Two;
gLikeThetaTwo = exp(-0.5*logLikeThetaTwo)*model.mixThetaTwo;
likeTheta = gLikeThetaOne + gLikeThetaTwo ...
    + (1-model.mixThetaOne-model.mixThetaTwo)/(2*pi);

thetaDiff = [0; thetaDiff; 0];
ratioOne = gLikeThetaOne./likeTheta;
ratioOne = [0; ratioOne; 0];
ratioTwo = gLikeThetaTwo./likeTheta;
ratioTwo = [0; ratioTwo; 0];
dL_dtheta = 1/model.sigma2One*...
    (-thetaDiff(1:end-1).*ratioOne(1:end-1) ...
     +thetaDiff(2:end).*ratioOne(2:end)) ...
    + 1/model.sigma2Two*...
    (-thetaDiff(1:end-1).*ratioTwo(1:end-1) ...
     +thetaDiff(2:end).*ratioTwo(2:end));


logLikeR1 = model.a*log(model.b) - gammaln(model.a) ...
    +(model.a-1)*log(model.r) - model.b*model.r + log(model.mixR);
logLikeR1minus = (model.a-1)*log(model.b) - gammaln(model.a-1) ...
    +(model.a-2)*log(model.r) - model.b*model.r + log(model.mixR);
logLikeR2 = log(model.b) - model.b*model.r + log(1-model.mixR);

logLikeR1(find(logLikeR1<-316))=-316;
logLikeR1minus(find(logLikeR1minus<-316))=-316;
logLikeR2(find(logLikeR2<-316))=-316;

ind = find(logLikeR1>logLikeR2);
logLikeR = zeros(size(logLikeR1));
logLikeR(ind) = log(1+exp(logLikeR2(ind))./exp(logLikeR1(ind)))+logLikeR1(ind);

ind2 = find(logLikeR1<=logLikeR2);
logLikeR(ind2) = log(1+exp(logLikeR1(ind2))./exp(logLikeR2(ind2)))+logLikeR2(ind2);

%ll = ll + sum(logLikeR);

%logLikeR = model.a*log(model.b) - gammaln(model.a) ...
%    +(model.a-1)*log(model.r) - model.b*model.r;
%logLikeRminus = (model.a-1)*log(model.b) - gammaln(model.a-1) ...
%    +(model.a-2)*log(model.r) - model.b*model.r;
%likeR = exp(logLikeR)*model.mixR;
likeR = exp(logLikeR); %likeR = likeR + (1-model.mixR)*model.b*exp(-model.b*model.r);

dL_dr = (model.b*(exp(logLikeR1minus) - exp(logLikeR1)) ...
         - model.b*exp(logLikeR2))./likeR;

X1 = model.X(1:end-1, :);
X2 = model.X(2:end, :);
diffX = X2 -X1;
r2 = model.r.*model.r;
ind = find(model.r);
a = diffX(:, 1)./model.r;
a2 = a.*a;
oneMinusa2 = 1-a2;
if any(oneMinusa2==0)
  oneMinusa2 = 1-a2+eps;
end
dr_dxi2 = diffX(:, 2)./model.r;
dr_dxi1 = diffX(:, 1)./model.r;
dr_dximin12 = -diffX(:, 2)./model.r;
dr_dximin11 = -diffX(:, 1)./model.r;

dTheta_dxi1 = (1./model.r-diffX(:, 1).*dr_dxi1./r2)./sqrt(oneMinusa2);
dTheta_dxi2 = (-diffX(:, 1).*dr_dxi2./r2)./sqrt(oneMinusa2);
dTheta_dximin11 = (-1./model.r-diffX(:, 1).*dr_dximin11./r2)./sqrt(oneMinusa2);
dTheta_dximin12 = (-diffX(:, 1).*dr_dximin12./r2)./sqrt(oneMinusa2);

gX(2:end, 2) = dL_dtheta.*dTheta_dxi2 + dL_dr.*dr_dxi2;
gX(2:end, 1) = dL_dtheta.*dTheta_dxi1 + dL_dr.*dr_dxi1;

gX(1:end-1, 2) = gX(1:end-1, 2) + dL_dtheta.*dTheta_dximin12 + dL_dr.*dr_dximin12;
gX(1:end-1, 1) = gX(1:end-1, 1) + dL_dtheta.*dTheta_dximin11 + dL_dr.*dr_dximin11;
