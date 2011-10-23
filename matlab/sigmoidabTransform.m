function y = sigmoidabTransform(x, transform, varargin)

% SIGMOIDABTRANSFORM Constrains a parameter to be between A and B
% by a scaled logistic sigmoid function.
%
% FORMAT
%
% DESC contains commands to constrain parameters to be between A
% and B via the sigmoid function, y=A+(B-A)/(1+exp(-x)).
%
% ARG x : input argument.
%
% ARG transform : type of transform, 'atox' maps a value into
% the transformed space (i.e. makes it between A and B). 'xtoa' 
% maps the parameter back from transformed space to the original
% space. 'gradfact' gives the factor needed to correct gradients
% with respect to the transformed parameter, that is, it gives
% the gradient dx/da where x is the transformed parameter and a
% is the untransformed parameter.
%
% ARG transformsettings (first element of varargin): vector [A B] 
% giving the minimum and maximum values A and B for the transformed 
% parameter. If not given, assume A=0 and B=1 as in the function 
% sigmoidTransform.
%
% OUTPUT y : return argument.
% 
% SEEALSO : negLogLogitTransform, expTransform
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007
%
% COPYRIGHT : Jaakko Peltonen, 2011

% OPTIMI


% Get the transformation settings from the varargin structure
if length(varargin)>0,
  transformsettings=varargin{1};
else
  transformsettings=[];
end;


identitymode=0;
if identitymode==1,
  % DEBUG ONLY, turns transformation to identity!
  fprintf(1,'Using identity mode for transform %s of parameter value %f\n',transform, x)
  switch transform
   case 'atox'
    y=x;
   case 'xtoa'
    y=x;
   case 'gradfact'
    y=0*x+1;
  end
else
  % normal mode, logistic sigmoid between A and B
  
  if isempty(transformsettings),
    warning(sprintf('sigmoidabTransform(%f,%s) called without transformation settings\n',x,transform));
    A=0; B=1;
  else
    A=transformsettings(1);
    B=transformsettings(2);
  end;
  
  limVal = 36;
  minval_sigmoid=A+(B-A)*eps;
  maxval_sigmoid=A+(B-A)*(1-eps);
  
  y = zeros(size(x));
  switch transform
   case 'atox'
    %x
    %A
    %B
    %fprintf(1,'transforming x %f, A %f, B %f\n',x,A,B);
    
    index = find(x<-limVal);
    y(index) = A+(B-A)*eps;
    x(index) = NaN;
    
    index = find(x<limVal);
    y(index) = A+(B-A)*sigmoid(x(index));
    x(index) = NaN;
    
    index = find(~isnan(x));
    y(index) = A+(B-A)*(1-eps);
    %y
    
   case 'xtoa'
    %  [x A B]
    %fprintf(1,'inverse-transforming x %f, A %f, B %f\n',x,A,B);
    index=find(x<=minval_sigmoid);
    y(index)=-limVal;
    x(index)=NaN;
    
    index=find(x<maxval_sigmoid);
    y(index) = invSigmoid((x(index)-A)/(B-A));
    x(index)=NaN;
    
    index=find(x>=maxval_sigmoid);
    y(index)=limVal;
    
   case 'gradfact'
    %fprintf(1,'gradfact-transforming x %f, A %f, B %f\n',x,A,B);
    %y = (B-A)*x.*(1-x);
    y = (x-A).*(B-x)/(B-A);
    %y=0*x+1; % DEBUG ONLY, turns factors off!
  end
end;

  
