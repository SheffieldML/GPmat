function pb = myProgressBar(arg1, arg2, counterPrintStep)

% MYPROGRESSBAR A text-based (no graphics) progress bar for long computations involving many iterations
% FORMAT
% DESC 
%
% COPYRIGHT : Andreas C. Damianou, 2013
%
% SEEALSO : myProgressBar2

% SHEFFIELDML

% Example: 
% pb = myProgressBar(length(testInd), length(testInd)/10); % This will print 10 times
% for i=1:length(testInd)
%    pb = myProgressBar(pb,i);
%    ...

if isstruct(arg1)
    pb = arg1;
    if nargin<2
        return
    end
    if mod(arg2,pb.counterPrintStep) == 0
        fprintf('.');
    end
else
    pb.totalIters = arg1;
    if nargin < 2
        pb.counterLength = 40;
    else
        pb.counterLength = arg2;
    end
    if nargin < 3
        pb.counterPrintStep = round(pb.totalIters/pb.counterLength);
    else
        pb.counterPrintStep = counterPrintStep;
    end
    
    
    for i=1:pb.counterLength
        fprintf(1,'.');
    end
    fprintf(1,'|\n');
end
