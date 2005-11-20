function y = gaussOverDiffCumGaussian(x, xp, order)

% GAUSSOVERDIFFCUMGAUSSIAN A gaussian in x over the difference between two cumulative Gaussians. 

% NDLUTIL

% Theoretically this is simply
% ngaussian(x)/(cumGaussian(x)-cumGaussian(xp)) but there are problems at
% extreme values. Order dictates whether ngaussian(x) or ngaussian(xp) is
% the numerator.

%.5*erfcx(-sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(x) ...
robustAdd = 1e-300;
fact = sqrt(2)/2;
xp2 = xp.*xp;
x2 = x.*x;
y = zeros(size(xp));
switch order
 case 1
  expRatio = exp(.5*(x2-xp2));
  index = find(x<=0);
  y(index) = 2./(erfcx(-fact*x(index)) ...
                 -expRatio(index).*erfcx(-fact*xp(index))+robustAdd); 
  x(index) = NaN;
  
  %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
  index = find(x>0);
  y(index) = 2./(expRatio(index).*erfcx(fact*xp(index))-erfcx(fact*x(index))+robustAdd); 
  y = y/sqrt(2*pi);
 case 2
  expRatio = exp(.5*(xp2-x2));
  index = find(x<=0);
  y(index) = 2./(expRatio(index).*erfcx(-fact*x(index)) ...
                 -erfcx(-fact*xp(index))+robustAdd); 
  x(index) = NaN;
  
  %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
  index = find(x>0);
  y(index) = 2./(erfcx(fact*xp(index))-expRatio(index).*erfcx(fact*x(index))+robustAdd); 
  y = y*1/sqrt(2*pi);
 otherwise
  error('Incorrect order')
end