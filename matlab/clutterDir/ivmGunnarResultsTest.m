function ivmGunnarResultsTest(dataSet);

% PPAGUNNARRESULTSTEST Helper script for collating results on Gunnar's benchmarks.

% PPA
beta = 0;
load(['ppa' dataSet 'Rbf']);
b = beta;
for dataNum = 1:10
  er(dataNum) = ppaGunnarTest(dataSet, dataNum, {'rbf', 'bias', 'white'}, ...
                                 kernParam, noiseParam, b);
end
save(dataSet, 'er');