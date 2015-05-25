function plot_expros(genes, exprs, exprsse, labels),

if nargin < 4,
  plot_exp(struct('m', getfields(exprs, genes), ...
		  'se', getfields(exprsse, genes)), genes);
else
  plot_exp(struct('m', getfields(exprs, genes), ...
		  'se', getfields(exprsse, genes)), labels);
end
