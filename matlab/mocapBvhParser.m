function [Y, connections] = mocapBvhParser(fileName)

% MOCAPBVHPARSER Parser of bvh files.

% MOCAP

lastToken = '';
jointNo = 1;

while(~feof(fid))
  nextLine = fgetl(fid);
  if ~isempty(nextLin)
    parts = stringSplit(nextLin, ' ');
    switch parts{1}
     case 'HIERARCHY'
      state = 'HIERARCHY';
     case 'MOTION'
      state = 'MOTION';
     otherwise
      switch state
       case 'HIERARCHY'
        switch parts{1}
         case {'ROOT', 'JOINT'}
          if strcmp(lastToken, '}')
            if length(joints) > 0
            end
          end
          
          end
         case 'MOTION'
        end 
      end
    end
    lastToken = parts{1};
     case 'ROOT'
     case 'OFFSET'
     case 'CHANNELS'
     case 'JOINT'
     case 'End'
     case '}'
     case '{'
     