function txt = escapeText(txt)
  
% ESCAPETEXT Add back slashes to escape existing backslashes in a text.

% SHEFFIELDML
  
txt = strrep(txt, '\', '\\');
txt = strrep(txt, '%', '%%');