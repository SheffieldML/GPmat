function y = gaussOverDiffCumGaussian(x, xp, order)

% GAUSSOVERDIFFCUMGAUSSIAN A Gaussian over difference of cumulative Gaussians.
% FORMAT
% DESC computes a Gaussian in x divided by the difference between
% two cumulative Gaussian distributions.
% ARG X1 : the argument of the first, positive, cumulative Gaussian.
% ARG X2 : the argument of the second, negative, cumulative Gaussian.
% ARG order : set to 1 or 2, determines whether X1 or X2 is used
% in the argument of the Gaussian term.
% RETURN y : returns y =
% ngaussian(X1)/(cumGaussian(X1)-cumGaussian(X2)) if order == 1 and
% ngaussian(X2)/(cumGaussian(X1)-cumGaussian(X2)) if order == 2
%
% Calculating this function naively causes problems at extreme values. 
%
% SEEALSO : lnCumGaussian, erfcx, lnDiffCumGaussian, cumGaussian
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% NDLUTIL


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
  error('Incorrect order, should be set to 1 or 2.')
end
