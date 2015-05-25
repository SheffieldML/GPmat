function answer = isoctave
  
% ISOCTAVE Returns true if the software running is Octave.
% FORMAT
% DESC tests if the software is octave or not.
% RETURN answer : true if the software is octave.
%
% COPYRIGHT : Neil D. Lawrence, 2008

% NDLUTIL
  
try 
  v = OCTAVE_VERSION;
  answer = true;
catch
  answer = false;
  return
end