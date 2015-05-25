function g = optimGradient(xvec, model, Y0, SoFC, sizePath, skel, varargin)

Xpath = reshape(xvec(1:2*sizePath), sizePath, 2);
Yfinal  = xvec(2*sizePath+1:2*sizePath + 62);

Yseq = zeros(sizePath,62);
Yseq(1,:)=Y0;
for i = 2:sizePath-1
    Yseq(i,:) = NaN;
end
Yseq(sizePath,:)=Yfinal;

g = -logLikeGradient(model, Xpath, Yseq);

finalXYZ   = convertxyz2Vect(SoFC);
initialXYZ = convertxyz2Vect(skel2xyz(skel,Yfinal));
for i=1:93

	if (finalXYZ(i)==0)
		tmp(i)        = 0;
	else
		tmp(i)        = 16*(finalXYZ(i) - initialXYZ(i));
	end
end
first      = computeGradient(skel,Yfinal);
g2  = 2*tmp*first;

[muFinal, covarSigma, factors] = gpPosteriorMeanCovar(model, Xpath(sizePath, :));

g2 = g2 - (pdinv(covarSigma)*(Yfinal - muFinal));

g(2*sizePath+1:2*sizePath + 62) = -g2;
