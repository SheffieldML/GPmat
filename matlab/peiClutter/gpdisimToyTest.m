% GPDISIMTOYTEST Script for toy test of GPDISIM model.

numGenes = 3;
numTimes = 10;
times = (1:numTimes)';

[genes, params] = gpdisimMakeToyData2(3, 10, 3);

for k = 1:length(genes),
  genevars = .01 * ones(size(genes{1}));
end

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
