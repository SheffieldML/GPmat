function [gradAlpha, gradOmega] = getLocalGradAlphaOmega(lfmKern)

% GETLOCALGRADALPHAOMEGA Gradients of parameters in alpha and omega
% FORMAT
% DESC computes gradients of alpha and omega, appearing in the LFM type of
% kernel, wrt parameters.
% ARG lfmKern : the kernel structure containing the parameters
% RETURN gradAlpha : gradients of parameters wrt alpha
% RETUNR gradOmega " gradinets of parameters wrt omega
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

m = lfmKern.mass;           % Par. 1
D = lfmKern.spring;         % Par. 2
C = lfmKern.damper;         % Par. 3
alpha = C/(2*m);
omega = sqrt(D/m-alpha^2);
% Derivatives of alpha and omega wrt parameters
gradAlpha = [-C/(2*m^2) 0 1/(2*m)];
gradOmega = [(C^2 - 2*m*D)/(4*m^3*omega) 1/(2*m*omega) -C/(4*m^2*omega)];
