function [xhist, Ehist, g_mean] = myhmc_bounded(findE, x, options, gradE, bounds, varargin),

L = options.iters;

Tau = options.tau;
epsilon = options.epsilon;

xhist = zeros(L, size(x, 2));
Ehist = zeros(L, 1);

if isfield(options, 'scales'),
  scales = options.scales;
else
  scales = diff(bounds) / max(diff(bounds));
end

M = diag(scales.^2);
M_invsqrt = diag(sqrt(1 ./ diag(M)));

g = feval(gradE, x, varargin{:}); % set gradient using initial x 
E = feval(findE, x, varargin{:}); % set objective function too 

g_mean = zeros(size(g));

for l = 1:L % loop L times 
  %p = mvnrnd(zeros(size(x')), M_inv, 1)' ; % initial momentum is Normal(0,1) 
  p = randn(size(x)) * M_invsqrt;
  H = p * M * p' / 2 + E ; % evaluate H(x,p) 
  xnew = x ; gnew = g ; 
  g_mean = g_mean + abs(g);
  try,
    for tau = 1:Tau % make Tau 'leapfrog' steps 
      p = p - epsilon * gnew / 2 ; % make half-step in p 
      xnew = xnew + epsilon * p * M ; % make step in x 
      while (any(xnew < bounds(1,:)) || any(xnew > bounds(2,:))),
	I1 = xnew < bounds(1, :);
	xnew(I1) = 2*bounds(1, I1) - xnew(I1);
	p(I1) = -p(I1);
	I2 = xnew > bounds(2, :);
	xnew(I2) = 2*bounds(2, I2) - xnew(I2);
	p(I2) = -p(I2);
      end
      gnew = feval(gradE, xnew, varargin{:}); % find new gradient 
      % if any(~isreal(gnew)),
      % 	fprintf('Unreal gradient!\n');
      % 	keyboard;
      % end
      p = p - epsilon * gnew / 2 ; % make half-step in p 
    end 
    Enew = feval(findE, xnew, varargin{:}) ; % find new value of H 
    Hnew = p * M * p' / 2 + Enew ; 
    dH = Hnew - H ; % Decide whether to accept 
    fprintf('Step %d, threshold: %.4f\n', l, exp(-dH));
    if ( dH < 0 ) accept = 1 ; 
    elseif ( rand() < exp(-dH) ) accept = 1 ; 
    else accept = 0 ; 
    end 
    if ( accept ) 
      fprintf('accepted.\n');
      if options.verbose > 1,
	disp(x);
      end
      g = gnew ; x = xnew ; E = Enew ; 
    end 
  catch,
    fprintf('Error in step %d\n', l);
  end
  xhist(l, :) = x;
  Ehist(l) = E;
end
g_mean = g_mean / L;
