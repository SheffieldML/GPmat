numGenes = 5;
numTimes = 12;
times = (1:numTimes)';

tinind = drosGetGeneinds(drosexp, drosTF.labels(1));
[origgenes, params] = gpdisimSimulateSystem(exp(reshape(drosexp.mean(tinind, :), ...
						  [12, 3])), ...
					5, 1:12);

for k = 1:length(origgenes),
  genes{k} = origgenes{k} ./ repmat(mean(origgenes{k}), [12, 1]);
  genes{k} = genes{k} + .1 * randn(size(genes{k}));
end
genevars = .01 * ones(size(genes{1}));

options = gpsimOptions;
model.type = 'cgpdisim';  % This new model type is a hack
for k=1:length(genes),
  model.comp{k} = gpdisimCreate(numGenes, 1, times, genes{k}, genevars, options);
end

options.fix(1).index = 2;
options.fix(1).value = expTransform(1, 'xtoa');
options.fix(2).index = 6;
options.fix(2).value = expTransform(1, 'xtoa');

model = modelOptimise(model, [], [], 1, 3000);
%model = gpdisimOptimise(model, 1, 3000);
