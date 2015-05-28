function    g = lfmGradientH(gamma1, gamma2, sigma2, gradThetaGamma, t1,t2,...
    mode, preComputeH, preComputeUpsilon)

% LFMGRADIENTH Gradient of the function h_i(z) with respect to some of the
%
%	Description:
%	hyperparameters of the kernel: m_k, C_k, D_k, m_r, C_r or D_r.
%
%	G = LFMGRADIENTH(GAMMA1, GAMMA2, SIGMA2, GRADTHETAGAMMA, T1, T2,
%	MODE) Computes the gradient of the function h_i(z) with respect to
%	some of the parameters of the system (mass, spring or damper).
%	 Returns:
%	  G - Gradient of the function with respect to the desired
%	   parameter.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  GRADTHETAGAMMA - Vector with the gradient of gamma1 and gamma2
%	   with respect to the desired parameter.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1)
%	  MODE - indicates in which way the vectors t1 and t2 must be
%	   transposed
%	
%	
%
%	See also
%	LFMKERNGRADIENT, LFMXLFMKERNGRADIENT, LFMGRADIENTUPSILON


%	Copyright (c) 2007, 2008 David Luengo
%	Copyright (c) 2008 Mauricio Alvarez



% Gradient evaluation

if gradThetaGamma(1) == 0
    if mode==1,
        % t1 is really t1 and t2 is really t2
        Tt1 = repmat(t1, 1, size(t2, 1));
    else
        % t1 is really t2 and t2 is really t1
        %Tt1 = repmat(t2', size(t1, 1), 1);
        Tt1 = repmat(t1', size(t2, 1), 1);
    end
    g = (gradThetaGamma(2)*Tt1.*exp(-gamma2*Tt1).* preComputeUpsilon ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
elseif gradThetaGamma(2) == 0
    if mode==1,
       % t1 is really t1 and t2 is really t2
       Tt1 = repmat(t1, 1, size(t2, 1));
         g = (lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2, t1,2) - exp(-gamma2*Tt1) ...
        .* lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2,zeros(size(t1)),4) ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
    else
        % t1 is really t2 and t2 is really t1
       Tt1 = repmat(t1', size(t2, 1),1);
       %Tt1 = repmat(t2', size(t1, 1), 1);
         g = (lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2, t1,1) - exp(-gamma2*Tt1) ...
        .* lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2,zeros(size(t1)),3) ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
    end
end
end




