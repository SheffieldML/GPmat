function y = gaussOverDiffCumGaussian1(x, xp, order)

% GAUSSOVERDIFFCUMGAUSSIAN1 A gaussian in x over the difference between two cumulative Gaussians. 

% IVM

% Theoretically this is simply
% ngaussian(x)/(cumGaussian(x)-cumGaussian(xp)) but there are problems at
% extreme values.

%.5*erfcx(-sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(x) ...
fact = sqrt(2)/2;
xp2 = xp.*xp;
x2 = x.*x;
y = zeros(size(xp));
switch order
  case 1
   expRatio = exp(.5*(x2-xp2));
   index = find(x<=0);
   y(index) = 2./(erfcx(-fact*x(index)) ...
                  -expRatio(index).*erfcx(-fact*xp(index))+eps); 
   x(index) = NaN;

   %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
   index = find(x>0);
   y(index) = 2./(expRatio(index).*erfcx(fact*xp(index))-erfcx(fact*x(index))+eps); 
   y = y/sqrt(2*pi);
 case 2
  expRatio = exp(.5*(xp2-x2));
  index = find(x<=0);
  y(index) = 2./(erfcx(-fact*expRatio.*x(index)) ...
                 -erfcx(-fact*xp(index))+eps); 
  x(index) = NaN;
  
  %.5*erfcx(sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(-x)=exp(.5*x*x)*(1-cumGaussian(x))
  index = find(x>0);
  y(index) = 2./(erfcx(fact*xp(index))-expRatio.*erfcx(fact*x(index))+eps); 
  y = y*1/sqrt(2*pi);
 otherwise
  error('Incorrect order')
end