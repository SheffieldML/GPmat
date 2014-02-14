function saveAllOpenFigures(path, prefix, closeAll)

% SAVEALLOPENFIGURES Save all open figures of the current MATLAB session. 
% DESC
% ARG path: the path where the figures are saved
% ARG prefix: prefix added to the saved files
% ARG closeAll: if set to true, the figures close after being saved
% 
% COPYRIGHT: Andreas C. Damianou, 2013
%
% SHEFFIELDML

% Don't forget to add the slash in the end, like: path = '../../TEMP/'

if ~exist('prefix')
    prefix=[];
end

figHandles = findobj('Type','figure');
 for i=1:length(figHandles)
    if ~exist(path)
        fprintf(1,'# saveAllOpenFigures: directory %s does not exist and is being created now...\n',path);
        
        % Create all the subdirectories from the root onwards, one by one
        seps = find(path=='/');
        subp=[];
        start = 1;
        for j=1:length(seps)
            last = seps(j);
            subp = [subp path(start:last)]
            start = seps(j)+1;
            com = ['mkdir ' subp];
            system(com);
        end
   
    end
	saveName = [path prefix num2str(i) '.png'];
	saveas(figHandles(i),saveName);
    fprintf(1,'* Saved %s\n',saveName);
 end

if exist('closeAll') && closeAll
    close all
end
