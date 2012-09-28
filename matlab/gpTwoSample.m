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
diffs(~diffs) = [];
typDiff = min(diffs);
invWidth = 1/(typDiff*typDiff);
%disp(1/invWidth)
options = gpOptions('ftc');
options.scale2var1 = true;
tfull = [t{1}; t{2}];
options.kern = kernCreate(tfull, options.kern);
options.kern.comp{1}.inverseWidth = invWidth;

optionsMulti = multimodelOptions('gp', 2, 'ftc');
optionsMulti.separate = [];
llComb = zeros(size(y{1}, 2), 1);
llInd = zeros(size(y{1}, 2), 1);
models = cell(size(y{1},2), 2);

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
    % Plot regression.
    %{
    xstar = linspace(t{1}(1), t{1}(end),100)';
    [mu, S] = gpPosteriorMeanVar(models{i,1}, xstar);
    % %       S = S - exp(2*loghypers(3,h)); % subtract noise variance
    f = [mu+2*sqrt(S); flipdim(mu-2*sqrt(S),1)];
    clf, hFill = fill([xstar; flipdim(xstar,1)], f, [0 0 1], 'FaceAlpha', .1, 'EdgeColor', [0 0 1], 'EdgeAlpha', 0.1);
    datamax = max(max(y{1,2})); datamin = min(min(y{1,2})); % Data min/max for plot limits.
    ylim([datamin datamax]), xlim([min(xstar) max(xstar)]),
    hold on, plot(xstar, mu,'b-','LineWidth',2)
    plot(t{1}, y{1}(:, i), '+b', t{2}, y{2}(:, i), 'or', 'MarkerSize', 8);
    set(get(get(hFill,'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');
    legend('GP','oxidised','non-oxidised')
    xlabel('time(months)','fontsize',12,'fontweight','bold'), ylabel('gene expression','fontsize',12,'fontweight','bold')
    %}
    
    % Fit independent Gaussian processes to the data.
    ymulti = cell(2, 1);
    ymulti{1} = y{1}(:, i);
    ymulti{2} = y{2}(:, i);
    optionsMulti.compOptions = options;
    multiModel = multimodelCreate(1, 1, t, ymulti, optionsMulti);
    % Make sure we scale using variance of combined data
    multiModel.comp{1}.scale = model.scale;
    multiModel.comp{1}.m = gpComputeM(multiModel.comp{1});
    multiModel.comp{2}.scale = model.scale;
    multiModel.comp{2}.m = gpComputeM(multiModel.comp{2});
    % Optimise the model
    multiModel.optimiser = 'scg';
    multiModel = modelOptimise(multiModel, [], [], display, maxIters, 0);
    llInd(i) = multimodelLogLikelihood(multiModel);
    models{i,2} = multiModel;
end
llratios = llInd - llComb;
end
