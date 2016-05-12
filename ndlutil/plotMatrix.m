function handle =  plotMatrix(A, ax, bracketStyle, type, options)

% PLOTMATRIX Fill a given axis with a matrix plot.
% FORMAT 
% DESC prints a matrix on the given axis.
% ARG A : matrix to print.
% ARG ax : axis handle.
% ARG bracketStyle : type of brackets to draw.
% ARG type : type of matrix to display.
% SEEALSO : printPlot
%
% COPYRIGHT : Neil D. Lawrence, 2010
  

% NDLUTIL

  if nargin < 5
      options = plotMatrixOptions;
      if nargin < 4
          if iscell(A)
              type = 'entries';
          else
              type = 'values';
          end
          if nargin < 3
              bracketStyle = 'square';
              if nargin < 2
                  ax = gca;
              end
          end
      end
  end
  if isfield(options, 'color') && ~isempty(options.color)
      lcolor = options.color;
  else
      lcolor = [0 0 0];
  end

  nrows = size(A, 1);
  ncols = size(A, 2);

  if strcmp('type', 'highlight')
      highlight = true;
      highlightRow = varargin{1};
      highlightCol = varargin{2};
      type = varargin{3};
      varargin = varargin{4:end};
      nargin = nargin - 3;
  end
      
  
  xlim = [-0.25 ncols+0.25]+0.5;
  xspan = xlim(2)-xlim(1);
  ylim = [-0.25 nrows+0.25]+0.5;
  yspan = ylim(2)-ylim(1);

  
  ca = gca;   
  axes(ax)
  cla
  switch type
    case 'image'
      handle =  image(A);

    case 'imagesc'
      handle =  imagesc(A, [min([min(min(A)) 0]) max(max(A))]);

    case 'values'
      for i = 1:nrows
          for j = 1:ncols
              handle(i, j) = text(j, i, num2str(A(i, j)), 'horizontalalignment', 'center');
          end
      end
    
    case 'entries'
      for i = 1:nrows
          for j = 1:ncols
              if ischar(A{i, j})
                  handle(i, j) = text(j, i, A{i, j}, ...
                                      'horizontalalignment', 'center');
              
              else              
                  handle(i, j) = text(j, i, ' ', ...
                                      'horizontalalignment', 'center');
              end
          end
      end
    case 'patch'
      for i = 1:nrows
          for j = 1:ncols
              handle(i, j) = patch([i - [0.5 0.5] i+[0.5 0.5]], ...
                    [j-0.5 j+0.5 j+0.5 j-0.5], (A(i, j))*[1 1 1]);
          end
      end
    case 'colorpatch'
      for i = 1:nrows
          for j = 1:ncols
              handle(i, j) = patch([i - [0.5 0.5] i+[0.5 0.5]], ...
                    [j-0.5 j+0.5 j+0.5 j-0.5], ...
                    [A(i, j, 1) A(i, j, 2) A(i, j, 3)]);
          end
      end
  end
  switch bracketStyle
    case 'boxes'
      for i = 0:nrows
          line([0 ncols] + 0.5, i + [0.5 0.5], 'color', lcolor);
      end
      for j = 0:ncols
          line(j + [0.5 0.5], [0 nrows]+0.5, 'color', lcolor);
      end
    case 'square'
      tickLength = 0.25;
      if isfield(options, 'bracketWidth') && ~isempty(options.bracketWidth)
          bracketWidth = options.bracketWidth;
      else
          bracketWidth = 3;
      end
      line([xlim(1)+tickLength xlim(1) xlim(1) xlim(1)+tickLength], ...
           [ylim(1) ylim(1) ylim(2) ylim(2)], ...
           'linewidth', bracketWidth, 'color', lcolor);
      line([xlim(2)-tickLength xlim(2) xlim(2) xlim(2)-tickLength], ...
           [ylim(1) ylim(1) ylim(2) ylim(2)], ...
           'linewidth', bracketWidth, 'color', lcolor);
    case 'none'
    otherwise
      error('Bracket type not specified')
  end
      
  if options.highlight.on
      if isfield(options.highlight, 'color') && ~isempty(options.highlight.color)
          hcolor = options.highlight.color;
      else
          hcolor = [0 0 0];
      end
      
      if isfield(options.highlight, 'width') && ~isempty(options.highlight.width)
          hwidth = options.highlight.width;
      else
          hwidth = 3;
      end


      
      
      hRow = options.highlight.row;
      hCol = options.highlight.col;
      if ischar(hRow) && hRow == ':'
          hRow = [1 nrows];
      end
      if ischar(hCol) && hCol == ':'
          hCol = [1 ncols];
      end
      if length(hRow) == 1
          hRow = [hRow hRow];
      end
      if length(hCol) == 1
          hCol = [hCol hCol];
      end
      hCol = sort(hCol);
      hRow = sort(hRow);
      line([hCol(1)-0.5 hCol(1)-0.5 ...
             hCol(2)+0.5 hCol(2)+0.5 hCol(1)-0.5], ...
            [hRow(1)-0.5 hRow(2)+0.5 ...
             hRow(2)+0.5 hRow(1)-0.5 hRow(1)-0.5], 'color', hcolor, 'linewidth', ...
           3)
  end

  if options.zoom.on
      
      zRow = options.zoom.row;
      zCol = options.zoom.col;
      if ischar(zRow) && zRow == ':'
          zRow = [1 nrows];
      end
      if ischar(zCol) && zCol == ':'
          zCol = [1 ncols];
      end
      if length(zRow) == 1
          zRow = [zRow zRow];
      end
      if length(zCol) == 1
          zCol = [zCol zCol];
      end
      zCol = sort(zCol);
      zRow = sort(zRow);
      xlim = [zCol(1)-0.5 zCol(2)+0.5];
      ylim = [zRow(1)-0.5 zRow(2)+0.5];
  end

  axis ij, axis equal, axis off
  set(gca, 'xlim', xlim);
  set(gca, 'ylim', ylim);

  if ~isempty(options.colormap)      
      colormap(options.colormap) 
  end
             
  axes(ca)
end
