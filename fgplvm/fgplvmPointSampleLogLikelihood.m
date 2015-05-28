function out = fgplvmPointSampleLogLikelihood(model,y,N,display)

% FGPLVMPOINTSAMPLELOGLIKELIHOOD
%
% COPYRIGHT : Carl Henrik Ek, 2008

% FGPLVM

if(nargin<4)
  display = true;
  if(nargin<3)
    N = 50;
    if(nargin<2)
      error('Too few Arguments');
    end
  end
end

if(model.q>3)
  error('Only Two and Three dimensional latent spaces supported');
end

% get limits
x_min = min(min(model.X));
x_max = max(max(model.X));

% get sample grid
switch model.q
 case 1
  G = linspace(x_min,x_max,N)';
 case 2
  [X Y] = meshgrid(linspace(x_min,x_max,N),linspace(x_min,x_max, ...
						    N));
  G = [reshape(X,prod(size(X)),1) reshape(Y,prod(size(Y)),1)];
 case 3
  [X Y Z] = meshgrid(linspace(x_min,x_max,N),linspace(x_min,x_max,N), ...
		     linspace(x_min,x_max,N));
  G = [reshape(X,prod(size(X)),1) reshape(Y,prod(size(Y)),1) reshape(Z,prod(size(Z)),1)];
end

l = zeros(size(G,1),1);
if(display)
  handle_waitbar = waitbar(0,'Computing point likelihood');
end
for(i = 1:1:size(G,1))
  l(i) = fgplvmPointLogLikelihood(model,G(i,:),y);
  if(display)
    waitbar(i/size(G,1));
  end
end
if(display)
  close(handle_waitbar);
end

switch model.q
 case 1
  out = l;clear l;
 case 2
  out = reshape(l,N,N);clear l;
 case 3
  out = reshape(l,N,N,N);clear l;
end

if(display)
  colormap gray;
  h = imagesc(out);colorbar;
end

return
