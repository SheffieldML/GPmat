function printLatexText(textString, fileName, directory)
  
% PRINTLATEXTEXT Print a text string to file for latex input.
% FORMAT 
% DESC prints a text string to file for latex input.
% ARG string : the string to print to file.
% ARG fileName : the base name of the file to use.
% ARG directory : the directory to place the eps files in.
% SEEALSO : printLatexPlot
%
% COPYRIGHT : Neil D. Lawrence, 2010
  
% NDLUTIL

  baseName = [directory filesep fileName];
  
  FID = fopen(baseName, 'w');    
  fprintf(FID, '%s', textString);
  fclose(FID);
    
end
