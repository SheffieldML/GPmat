function pb = myProgressBar2(pb, arg1, customMsg)

% MYPROGRESSBAR2 A graphics-based (progress bar for long computations involving many iterations
% FORMAT
% DESC 
%
% COPYRIGHT : Andreas C. Damianou, 2013
%
% SEEALSO : myProgressBar

% SHEFFIELDML

if isempty(pb) % first call
    totalIt = arg1;
    pb.curIt = 1;
    pb.totalIt = totalIt;
    pb.start = tic;
    pb.h = waitbar(0,'Initializing waitbar...');
else
    if nargin > 1 && ~isempty(arg1)
        incrIters = arg1;
    else
        incrIters = 1;
    end
    if nargin < 3
        customMsg = [];
    end
    pb.curIt = pb.curIt + incrIters;
    pb.timeSoFar = toc(pb.start);
    pb.totalTime = pb.totalIt * pb.timeSoFar / pb.curIt;
    pb.timeRemaining = pb.totalTime - pb.timeSoFar;
    waitbar(pb.curIt / pb.totalIt, pb.h, [secs2hms(pb.timeRemaining) ' remaining' customMsg]); 
    %fprintf('# Time remaining (min): \b %.2f', pb.timeRemaining/60)
end
