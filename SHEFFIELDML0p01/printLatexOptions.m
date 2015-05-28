function options = printLatexOptions

% PRINTLATEXOPTIONS Options for printing a plot to LaTeX.
%
%	Description:
%
%	OPTIONS = PRINTLATEXOPTIONS
%	 Returns:
%	  OPTIONS - options specify whether to create an eps and keep aspect
%	   ratio. SEEALSO : printLatexPlot


%	Copyright (c) 2010 Neil D. Lawrence

  
  options.pdf = true;
  options.eps = false;
  options.maintainAspect = true;
  options.height = 0;
  options.backgroundFile = '';
  options.axisColor = [];
end