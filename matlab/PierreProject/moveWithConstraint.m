%1 18 35 62

model.dynamics.kern.comp{2}.variance=0.1;

start = 18;
XinS = 18;
XinF = 35;
iter = 100;

SoFC = zeros(31,3);
%Feet on the ground

SoFC(28,:) = [6,14,0];
SoFC(21,:) = [14,8,0];


Xstart = model.X(start,:);
Y0 = model.y(start,:);


sizePath = XinF - XinS + 1;
Xin = model.X(XinS:XinF,:);

xvec = Xin(:)';
xvec(2*sizePath+1:2*sizePath+62) = model.y(XinF, :);

Xout = findPathOpt(model, xvec, Y0, SoFC, sizePath, skel, 'true', iter);

Xoutr = reshape(Xout(1:2*sizePath), sizePath, 2);
Yfinal = Xout(2*sizePath+1:2*sizePath+62);

pause;
figure(2)
YPath = fgplvmPosteriorMeanVar(model, Xoutr);
skelPlayData(skel, YPath,1/30);
plot3(SoFC(:,1), SoFC(:,3), SoFC(:,2), 'ko');

XYZFinal = skel2xyz(skel,Yfinal)
XYZPathFinal = skel2xyz(skel,YPath(sizePath,:))


%figure(3);
hold on;
skelVisualise(Yfinal, skel);
hold on;
plot3(XYZPathFinal(:,1), XYZPathFinal(:,3), XYZPathFinal(:,2), 'bo');
hold on;
plot3(XYZFinal(:,1), XYZFinal(:,3), XYZFinal(:,2), 'r+');
hold off;