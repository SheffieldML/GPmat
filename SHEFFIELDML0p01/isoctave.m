function answer = isoctave

% ISOCTAVE Returns true if the software running is Octave.
%
%	Description:
%
%	ANSWER = ISOCTAVE tests if the software is octave or not.
%	 Returns:
%	  ANSWER - true if the software is octave.


%	Copyright (c) 2008 Neil D. Lawrence

  
try 
  v = OCTAVE_VERSION;
  answer = true;
catch
  answer = false;
  return
end