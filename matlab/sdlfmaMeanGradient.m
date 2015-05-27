function [gradAlpha, gradOmega] = sdlfmaMeanGradient(sdlfmKern, t , option)

% SDLFMAMEANGRADIENT Gradients wrt parameters of the accel. mean SDLFM.
% Computes the gradients of the terms $r_d$ and $b_d$ that appear in the 
% mean function with respect to $\alpha_d$ and $\omega_d$. 
% FORMAT
% DESC computes the gradients of the terms $r_d$ that appear in the 
% mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN gradAlpha : gradient of $r_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $b_d$ wrt $\omega_d$
%
% FORMAT
% DESC computes the gradients of the terms $r_d$ and $b_d$ that appear in 
% the mean function with respect to $\alpha_d$ and $\omega_d$. 
% ARG sdlfmKern : switching dynamical LFM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% ARG option : indicates which to which term of the mean should be should
% the derivatives wrt $\alpha_d$ and  $\omega_d$ should be computed. Option
% 'Pos' computes the respective derivatives of the term $r_d$ and option 
% 'Vel' computes the respective derivatives of $b_d$.
% RETURN gradAlpha : gradient of $r_d$ or $b_d$ wrt $\alpha_d$
% RETURN gradOmega : gradient of $r_d$ or $b_d$ wrt $\omega_d$
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 3
    option = 'Pos';
end
    
alpha = sdlfmKern.damper/(2*sdlfmKern.mass);
omega = sqrt(sdlfmKern.spring/sdlfmKern.mass-alpha^2);
freq = omega*t;
hd = sdlfmvMeanCompute(sdlfmKern, t, 'Vel');

switch option
    case 'Pos'
        gd = sdlfmvMeanCompute(sdlfmKern, t);
        rd = sdlfmaMeanCompute(sdlfmKern, t);
        gradAlpha = -t.*rd - 2*alpha*hd - gd;
        [void, gradOmegagd] = sdlfmvMeanGradient(sdlfmKern, t);
        gradOmega = -alpha*gradOmegagd - omega*exp(-alpha*t).*cos(freq)*(1- alpha^2/omega^2) ...
            + (alpha^2/omega + omega)*exp(-alpha*t).*(-cos(freq) + omega*t.*sin(freq));
              
    case 'Vel'
        bd = sdlfmaMeanCompute(sdlfmKern, t, 'Vel');        
        gradAlpha = -t.*bd - 2*hd;
        gradOmega = exp(-alpha*t).*((-1- alpha^2/omega^2)*sin(freq) + ...
            (alpha^2/omega - omega)*t.*cos(freq) + 2*alpha*t.*sin(freq));
    otherwise
     error('No recognized option')   
end
