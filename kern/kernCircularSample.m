function includeText = kernCircularSample(mu, K, numSamps, numTheta, diagInfo)

% KERNCIRCULARSAMPLE Sample from covariance in a circular way alla Hennig.
%
% FORMAT
% DESC samples from GP along elipses of equiprobability
% ARG mu : mean of GP.
% ARG K : covariance of GP.
% ARG numSamps : number of samples from the Gaussian process.
% ARG diagInfo : information for the diagram.
% RETURN includeText : the text used to include in latex.
%
% SEEALSO : kernCreate
%
% COPYRIGHT: Neil D. Lawrence, 2013

% GPMAT

tau = 2*pi;

N = size(K, 1);
numSamps = length(diagInfo.lh);

thetaVector = linspace(0, tau, numTheta+1);
thetaVector = thetaVector(1:end-1);

R1 = randn(N, numSamps);
U1 = R1*diag(1./sqrt(sum(R1.*R1, 1)));
R2 = randn(N, numSamps);
R2 = R2 - U1*diag(sum(R2.*U1, 1));
R2 = R2*diag(sqrt(sum(R1.*R1, 1))./sqrt(sum(R2.*R2, 1)));

L = jitChol(K)';

includeText = '';
count = 0;
for theta = thetaVector
  count = count + 1;
  xc = cos(theta);
  yc = sin(theta);
  % generate 2d basis in t-d space
  coord = xc*R1 + yc*R2;
  y = L*coord;
  for i = 1:numSamps
    set(diagInfo.lh(i), 'ydata', y(:, i));
  end
  fileName = [diagInfo.fileBase num2str(count)];
  printLatexPlot(fileName, diagInfo.dirName, diagInfo.colWidth);
  if count > 1
    includeText = [includeText '\newframe'];
  end
  includeText = [includeText '\inputdiagram{' diagInfo.dirName fileName '}'];
end
printLatexText(includeText, [diagInfo.fileBase 'IncludeText.tex'], diagInfo.dirName);
