function term  = sdlfmjMeanCompute(sdlfmKern, t , option)

% SDLFMJMEANCOMPUTE Jolt mean for the switching dynamical LFM model.
%
%	Description:
%	Computes the terms $r_d$ and $q_d$ that appear in the mean function
%	associated with the switching dynamical LFM model. If the mean function
%	is mu(t), then
%	
%	mu(t) = r_d(t)y_d(t_0) + q_d(t)\dot{y}_d(t_0),
%	
%	where $y_d(t_0)$ is the initial condition associated to the position and
%	$\dot{y}_d(t_0)$ is the initial condition associated to the velocity.
%	
%
%	TERM = SDLFMJMEANCOMPUTE(SDLFMKERN, T) 
%	 Returns:
%	  TERM - the value of $r_d$.
%	 Arguments:
%	  SDLFMKERN - switching dynamical LFM kernel structure with the
%	   parameters.
%	  T - input times for which the mean is to be computed.
%
%	TERM = SDLFMJMEANCOMPUTE(SDLFMKERN, T, OPTION)  Computes the terms
%	that appear in the mean function associated with the switching
%	dynamical LFM model.
%	 Returns:
%	  TERM - the value of $r_d$ or $q_d$ depending of option.
%	 Arguments:
%	  SDLFMKERN - switching dynamical LFM kernel structure with the
%	   parameters.
%	  T - input times for which the mean is to be computed.
%	  OPTION - indicates which term of the mean should be computed.
%	   Option 'Pos' computes the term $r_d$ and option 'Vel' computes
%	   $q_d$ that accompanies the initial condition of the velocity.


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin < 3
    option = 'Pos';
end
    
alpha = sdlfmKern.damper/(2*sdlfmKern.mass);
omega = sqrt(sdlfmKern.spring/sdlfmKern.mass-alpha^2);
freq = omega*t;

switch option
    case 'Pos'       
        term = (alpha^2/omega + omega)*exp(-alpha*t).*...
            ((omega^2 - alpha^2)*sin(freq) + 2*alpha*omega*cos(freq));
    case 'Vel'
        term = exp(-alpha*t).*((3*alpha*omega - alpha^3/omega)*sin(freq)...
            +(3*alpha^2- omega^2)*cos(freq));         
    otherwise
     error('No recognized option')   
end
