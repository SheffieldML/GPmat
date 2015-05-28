% DEMCMU35NEARESTNEIGHBOUR Recreate the Nearest Neighbour result from Taylor et al.

datasetName = 'cmu35Taylor';
[Y, lbls, Ytest, lblsTest] = lvmLoadData(datasetName);
load cmu35TaylorScaleBias
bias = [bias bias];
scale = [scale scale];
Y = Y - repmat(bias, size(Y, 1), 1);
Ytest = Ytest - repmat(bias, size(Ytest, 1), 1);
Y = Y.*repmat(scale, size(Y, 1), 1);
Ytest = Ytest.*repmat(scale, size(Ytest, 1), 1);

% Indices associated with right leg.
legInd = [17:23 76:82];
leglessInd = [1:16 24:75 83:118];
YtrainLeg = Y(:, leglessInd);
jointAnglesInd = [19:size(YtrainLeg, 2)];

startInd = 63;
YtestLeg = Ytest(startInd:end, leglessInd);

dists = dist2(YtrainLeg(:, jointAnglesInd), YtestLeg(:, jointAnglesInd));
[void, bestIndLeg] = min(dists);
Ydiff = (Ytest(startInd:end, legInd) - Y(bestIndLeg, legInd));
errLeg = sum(sum(Ydiff.*Ydiff))/length(legInd);


figure(1)
clf
subplot(1, 2, 1);
plot(45:size(Ytest, 1), Ytest(45:end, 8), '-')
hold on
plot(45:size(Ytest, 1), [Ytest(45:startInd-1, 8); Y(bestIndLeg, 8)], '--')
subplot(1, 2, 2);
plot(45:size(Ytest, 1), Ytest(45:end, 9), '-')
hold on
plot(45:size(Ytest, 1), [Ytest(45:startInd-1, 9); Y(bestIndLeg, 9)], '--')


% Indices associated with upper body.
bodyInd = [30:59 89:118];
bodylessInd = [1:29 60:88];
YtrainBody = Y(:, bodylessInd);
jointAnglesInd = [19:size(YtrainBody, 2)];
startInd = 63;
YtestBody = Ytest(startInd:end, bodylessInd);

dists = dist2(YtrainBody(:, jointAnglesInd), YtestBody(:, jointAnglesInd));
[void, bestIndBody] = min(dists);
Ydiff = (Ytest(startInd:end, bodyInd) - Y(bestIndBody, bodyInd));
errBody = sum(sum(Ydiff.*Ydiff))/length(bodyInd);

figure(1)
clf
subplot(1, 2, 1);
plot(45:size(Ytest, 1), Ytest(45:end, 49), '-')
hold on
plot(45:size(Ytest, 1), [Ytest(45:startInd-1, 49); Y(bestIndBody, 8)], '--')
subplot(1, 2, 2);
plot(45:size(Ytest, 1), Ytest(45:end, 50), '-')
hold on
plot(45:size(Ytest, 1), [Ytest(45:startInd-1, 50); Y(bestIndBody, 9)], '--')


