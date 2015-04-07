function r = plot_tf_exps(exp_struct, expse_struct),

TFs = {'tin', 'bin', 'twi', 'bap', 'mef2'};
TFFBs = {'FBgn0004862', 'FBgn0045759', 'FBgn0011656', ...
	 'FBgn0004110', 'FBgn0003900'};

if nargin > 1,
  r.m = zeros(5, 36);
  r.se = zeros(5, 36);

  for k=1:5,
    r.m(k, :) = exp_struct.(TFFBs{k});
    r.se(k, :) = expse_struct.(TFFBs{k});
  end
else
  r = zeros(5, 36);

  for k=1:5,
    r(k, :) = exp_struct.(TFFBs{k});
  end
end

plot_exp(r, TFs);
