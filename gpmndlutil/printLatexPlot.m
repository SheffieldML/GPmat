function printLatexPlot(fileName, directory, width, options)
  
% PRINTLATEXPLOT Print a plot to LaTeX.
% FORMAT 
% DESC prints a plot to the specified file name and directory.
% ARG fileName : the base name of the file to use.
% ARG directory : the directory to place the eps files in.
% ARG width : the width in cm of the plot (default 6).
% ARG options : options specify whether to create an eps and keep aspect ratio. 
% SEEALSO : printPlot
%
% COPYRIGHT : Neil D. Lawrence, 2010
  
% NDLUTIL

  if nargin<4
    options = printLatexOptions;
    if nargin < 3
      width = 6;
    end
  end

  pos = get(gcf, 'paperposition');
  origpos = pos;
  if options.maintainAspect == true
    height = pos(4)/pos(3)*width;
  else 
    height = options.height;
  end
  global printDiagram

  if isempty(printDiagram)
    printTrue = true; % User hasn't been explicit, assume print.
  else
    printTrue = printDiagram; % Use user preference.
  end
  if nargin < 2 
    directory = '.';
  end
  if nargin < 3
    png = false;
  else
    png = true;
  end
  if isoctave
    pos(3) = width/2.54;
    pos(4) = height/2.54;
    set(gcf, 'paperposition', pos);
    baseName = [pwd filesep fileName];
    newName = [pwd filesep directory filesep fileName];
    fprintf('Printing LaTeX plot ...\n');
    if ~options.tikz
      print('-dtex', [fileName '.tex']);
      if options.pdf
        system(['epstopdf ' baseName '.eps --outfile ' baseName '.pdf']);
        system(['mv ' baseName '.pdf ' newName '.pdf']);
        if ~options.eps
          system(['rm ' baseName '.eps']);
        else
          system(['mv ' baseName '.eps ' newName '.eps']);
        end
      end
    
      % Use another background file. 
      if length(options.backgroundFile) == 0
        backgroundFile = fileName;
      else
        if options.pdf
          system(['rm ' baseName '.pdf'])
        end
        if options.eps
          system(['rm ' baseName '.eps'])
        end
        backgroundFile = options.backgroundFile;
      end
    
      FID = fopen([baseName '.tex'], 'r');
      fileStr = '';
      lin = getline(FID);
      i = 0;
      while(lin>0)
        i = i + 1;
        lin = regexprep(lin, ['\\includegraphics{' fileName '}'], ...
                        ['\includegraphics{' directory filesep backgroundFile '}']);
        fileStr = [fileStr escapeText(lin) '\n'];
        %disp(fileStr)
        lin = getline(FID);
      end
      fclose(FID);
      system(['rm -f ' fileName '.tex']);
      FID = fopen([newName '.tex'], 'w');
      fprintf(FID, fileStr);
      fclose(FID);
    else
      print('-dtikz', [newName '.tex']);
    end
  else
    matfig2pgf('filename', [directory filesep backgroundFile '.tex'], ...
               'figwidth', width)
  end
end
