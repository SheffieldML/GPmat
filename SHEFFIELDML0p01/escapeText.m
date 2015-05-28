function txt = escapeText(txt)

% ESCAPETEXT Add back slashes to escape existing backslashes in a text.
%
%	Description:
%	txt = escapeText(txt)
%
  
txt = strrep(txt, '\', '\\');
txt = strrep(txt, '%', '%%');