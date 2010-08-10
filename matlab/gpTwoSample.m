function [llratios, models] = gpTwoSample(t, y)
  
% GPTWOSAMPLE Do Oliver Stegles simple two sample test on a data set.
% FORMAT
% DESC computes a two sample test using independent or combined GPs to
% determine whether two gene expression profiles were more likely to be
% independently generated or generated together.
% ARG t : cell array containing times of the two gene expression
% profiles.
% ARG y : cell array containing log expression levels of the gene
% expression. In each component of the cell array different genes are in
% different columns and different time points in different rows.
% RETURN llratios : the log likelihood ratio between the combined and
% independent models.
%
% COPYRIGHT : Neil D. Lawrence, 2010
%
% SEEALSO : multimodelCreate, gpCreate
  
% GP
  maxIters = 500;
  display = 0;
  tTemp = sort(t{1});
  diffs = tTemp(2:end) - tTemp(1:end-1);
  diffs(find(~diffs)) = [];
  typDiff = min(diffs);
  invWidth = 1/(typDiff*typDiff);
  options = gpOptions('ftc');
  options.scale2var1 = true;
  tfull = [t{1}; t{2}];
  options.kern = kernCreate(tfull, options.kern);
  options.kern.comp{1}.inverseWidth = invWidth;

  optionsMulti = multimodelOptions('gp', 2, 'ftc');
  optionsMulti.separate = [];
  llComb = zeros(size(y{1}, 2), 1);
  llInd = zeros(size(y{1}, 2), 1);

  for i = 1:size(y{1}, 2)
    fprintf('Start gene %d\n', i)

    % Check for missing data.
    if any(isnan(y{1}(:, i))) || any(isnan(y{2}(:, i)))
      options.isMissingData = true;
      options.isSpherical = false;
    else
      options.isMissingData = false;
      options.isSpherical = true;
    end
    % Fit a combined Gaussian process to the data.
    model = gpCreate(1, 1, tfull, [y{1}(:, i); y{2}(:, i)], options);
    
    model = gpOptimise(model, display, maxIters, 0);
    llComb(i) = gpLogLikelihood(model);
    models{i, 1} = model;

    % Fit independent Gaussian processes to the data.
    ymulti = cell(2, 1);
    ymulti{1} = y{1}(:, i);
    ymulti{2} = y{2}(:, i);
    optionsMulti.compOptions = options;
    multiModel = multimodelCreate(1, 1, t, ymulti, optionsMulti);

    % Make sure we scale using variance of combined data
    multiModel.comp{1}.scale = model.scale;
    multiModel.comp{2}.scale = model.scale;

    % Optimise the model
    multiModel.optimiser = 'scg';
    multiModel = modelOptimise(multiModel, [], [], display, maxIters, 0);
    llInd(i) = multimodelLogLikelihood(multiModel);
    models{i, 2} = multiModel;
  end
  llratios = llInd - llComb;
end
