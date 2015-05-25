function f = optimSequence(xvec, model, Y0, SoFC, sizePath, skel, varargin)

Xpath = reshape(xvec(1:2*sizePath), sizePath, 2);
Yfinal  = xvec(2*sizePath+1:2*sizePath + 62);

Yseq = zeros(sizePath,62);
Yseq(1,:)=Y0;
for i = 2:sizePath-1
    Yseq(i,:) = NaN;
end
Yseq(sizePath,:)=Yfinal;

A = -logLikelihood(model, Xpath, Yseq);
 
finalXYZ   = convertxyz2Vect(SoFC);
initialXYZ = convertxyz2Vect(skel2xyz(skel,Yfinal));

for i=1:93

	if (finalXYZ(i)== 0)
		tmp(i)        = 0;
	else
		tmp(i)        = 4*(finalXYZ(i) - initialXYZ(i));
	end
end



[muFinal, covarSigma, factors] = gpPosteriorMeanCovar(model, Xpath(sizePath, :));
diff = (Yfinal - muFinal)';
B= 0.5*(diff'*pdinv(covarSigma)*diff);

for i=1:size(tmp,2)
    B = B + tmp(1,i)*tmp(1,i);
end
f =  A + B;
%animOptim(model, xvec, sizePath);
