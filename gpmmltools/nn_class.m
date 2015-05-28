function [class, d_class] = nn_class(X,X_test,k,dist_type) 

% NN_CLASS
%
%
%

% SHEFFIELDML

if(nargin<4)
  dist_type = 'euclidean';
  if(nargin<3)
    k = 1;
    if(nargin<2)
      error('To Few Arguments');
    end
  end
end

dist_func = str2func(strcat('dist_',dist_type));
d = dist_func(X,X_test);

class = zeros(size(X_test,1),k);
d_class = zeros(size(class));
for(i = 1:1:size(X_test,1))
  [tmp ind] = sort(d(:,i),'ascend');
  class(i,:) = ind(1:1:k);
  d_class(i,:) = tmp(1:1:k);
end

return

%------------------------------------------%

function D = dist_euclidean( X, Y )
if( ~isa(X,'double') || ~isa(Y,'double'))
  error( 'Inputs must be of type double'); 
end;
m = size(X,1); n = size(Y,1);  
Yt = Y';  
XX = sum(X.*X,2);        
YY = sum(Yt.*Yt,1);      
D = XX(:,ones(1,n)) + YY(ones(1,m),:) - 2*X*Yt;

return

function D = dist_emd( X, Y )
[m p] = size(X);  [n p] = size(Y);  

Xcdf = cumsum(X,2);
Ycdf = cumsum(Y,2);

m_ones = ones(1,m); D = zeros(m,n);
for i=1:n
  ycdf = Ycdf(i,:); 
  ycdf_rep = ycdf( m_ones, : );
  D(:,i) = sum(abs(Xcdf - ycdf_rep),2);
end

return

function D = dist_chisquared( X, Y )
[m p] = size(X);  [n p] = size(Y);

m_ones = ones(1,m); D = zeros(m,n); 
for i=1:n  
  yi = Y(i,:);  yi_rep = yi( m_ones, : );
  s = yi_rep + X;    d = yi_rep - X;
  D(:,i) = sum( d.^2 ./ (s+eps), 2 );
end
D = D/2;

return

function D = dist_L1( X, Y )
[m p] = size(X);  [n p] = size(Y);

m_ones = ones(1,m); D = zeros(m,n); 
for i=1:n  
  yi = Y(i,:);  yi = yi( m_ones, : );
  D(:,i) = sum( abs( X-yi),2 );
end

return

