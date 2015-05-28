function gX =  whiteXwhiteKernGradX(whiteKern1, whiteKern2, X, X2) 

% WHITEXWHITEKERNGRADX
  
% KERN

if nargin < 3,    
    X2 = X;
else
    U = X;
    X = X2;
    X2 = U;
end
    

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
