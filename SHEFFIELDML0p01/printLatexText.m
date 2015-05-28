function printLatexText(textString, fileName, directory)

% PRINTLATEXTEXT Print a text string to file for latex input.
%
%	Description:
%
%	PRINTLATEXTEXT(STRING, FILENAME, DIRECTORY) prints a text string to
%	file for latex input.
%	 Arguments:
%	  STRING - the string to print to file.
%	  FILENAME - the base name of the file to use.
%	  DIRECTORY - the directory to place the eps files in. SEEALSO :
%	   printLatexPlot


%	Copyright (c) 2010 Neil D. Lawrence


  baseName = [directory filesep fileName];
  
  FID = fopen(baseName, 'w');    
  fprintf(FID, '%s', textString);
  fclose(FID);
    
end
