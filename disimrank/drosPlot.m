function drosPlot(model, totgenes, thisgene)

% DROSPLOT Plot a GPDISIM model/models
% FORMAT
% DESC Plot a GPDISIM model/models
% ARG model : The model to plot
% ARG totgenes : Total number of genes if plotting several single
% target models in a single figure (optional)
% ARG thisgene : Index of this gene in the above case (optional)
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

FONTSIZE = 8;
LINEWIDTH = 1;
MARKERSIZE = 6;

if strcmp(model.type, 'cgpdisim'),
  numGenes = model.comp{1}.numGenes + 1;
else
  numGenes = model.comp{1}.numGenes;
end

selectgenes = 1:numGenes;

numPlots = length(selectgenes);

tf = model.comp{1}.annotation.tf;
genes = model.comp{1}.annotation.genes;
targets = genes(2:end);

genenames = genes;
times = model.comp{1}.t;

y = {};
yvar = {};
for k=1:length(model.comp),
  y{k} = reshape(model.comp{k}.y, [length(times), length(genes)]);
  yvar{k} = reshape(model.comp{k}.yvar, [length(times), length(genes)]);
end

for j = 1:length(model.comp)

  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).
  predt = [1:0.1:12 model.comp{j}.t']';

  if model.comp{j}.includeNoise,
    simMultiKern = model.comp{j}.kern.comp{1};
  else
    simMultiKern = model.comp{j}.kern;
  end
  if strcmp(model.type, 'cgpdisim'),
    proteinKern = kernCreate(model.comp{j}.t, 'sim');
    proteinKern.inverseWidth = simMultiKern.comp{1}.inverseWidth;
    proteinKern.decay = model.comp{j}.delta;
    proteinKern.variance = simMultiKern.comp{2}.di_variance;  
    % simXrbf requires an RBF with unit variance
    inputKern = kernCreate(model.comp{j}.t, 'rbf');
    inputKern.inverseWidth = simMultiKern.comp{1}.inverseWidth;
    inputKern.variance = 1;
    K = simXrbfKernCompute(proteinKern, inputKern, ...
			   predt, model.comp{j}.t);
    K = K * simMultiKern.comp{1}.variance;
    warning('off', 'KERN:simRBFVariance');
    for i=2:simMultiKern.numBlocks
      blockK =  disimXsimKernCompute(simMultiKern.comp{i}, proteinKern, ...
				     model.comp{j}.t, predt);
      K = [K blockK'];    
    end
    warning('on', 'KERN:simRBFVariance');
  
    ymean = reshape(ones(length(times),1)*[0, model.comp{j}.mu], ...
		    size(model.comp{j}.y));
  else
    proteinKern = kernCreate(model.comp{1}.t, 'rbf'); 
    proteinKern.inverseWidth = ...
	simMultiKern.comp{1}.inverseWidth;
    K = [];
  
    for i=1:simMultiKern.numBlocks
      K = [K; simXrbfKernCompute(simMultiKern.comp{i}, proteinKern, ...
				 model.comp{j}.t, predt)];
    end
    
    K = K';
    ymean = reshape(ones(length(times),1)*model.comp{j}.mu, ...
		    size(model.comp{j}.y));
  end
  
  predF = K*model.comp{j}.invK*(model.comp{j}.y-ymean);
  if strcmp(model.type, 'cgpdisim'),
    varF = simMultiKern.comp{1}.variance * ...
	   kernDiagCompute(proteinKern, predt) - sum(K'.* ...
						     (model.comp{j}.invK* ...
						      K'), 1)';
  else
    varF = kernDiagCompute(proteinKern, predt) - sum(K'.* ...
						     (model.comp{j}.invK* ...
						      K'), 1)';
  end

  % Predicted Gene Expressions
  if strcmp(model.type, 'cgpdisim'),
    Kxx = multiKernCompute(simMultiKern, predt, model.comp{j}.t);
    varX = real(kernDiagCompute(simMultiKern, predt) - ...
		sum(Kxx'.* (model.comp{j}.invK*Kxx'), 1)');
    meanPredX = reshape(ones(length(predt),1)*([0 model.comp{j}.B./ ...
		    model.comp{j}.D]), length(predt)*(numGenes), 1);
  else
    Kxx = multiKernCompute(simMultiKern, predt, model.comp{j}.t);
    varX = real(kernDiagCompute(simMultiKern, predt) - ...
		sum(Kxx'.* (model.comp{j}.invK*Kxx'), 1)');
    meanPredX = reshape(ones(length(predt),1)*(model.comp{j}.B./ ...
		    model.comp{j}.D), length(predt)*(numGenes), 1);
  end

  predX = meanPredX + real(Kxx*model.comp{j}.invK*(model.comp{j}.y-ymean));

  % Take out predictions at data points.
  % Use them to get the scale for the other data.
  numData = length(predX)/(numGenes);
  predExprs = reshape(predX, numData, numGenes);
  meanExprs = ones(numData, 1)*mean(predExprs);
  scaleExprs = ones(numData, 1)*(sqrt(var(predExprs))./sqrt(var(y{j})));

  predExprs(end-length(model.comp{j}.t)+1:end,:) = [];
  varExprs = reshape(varX, numData, numGenes);
  varExprs = varExprs./scaleExprs./scaleExprs;
  varExprs(end-length(model.comp{j}.t)+1:end,:) = []; 
  predF(end-length(model.comp{j}.t)+1:end) = [];
  varF(end-length(model.comp{j}.t)+1:end) = [];
  predt(end-length(model.comp{j}.t)+1:end) = [];
  
  yscale = ones(length(times), 1)*sqrt(var(y{j}));
  ymean = ones(length(times),1)*mean(y{j});
  ynormal = ymean + (y{j}-ymean)./yscale;
  yvarNormal = yvar{j}./yscale./yscale;

  scalePred = sqrt(var(predExprs));
  
  if nargin < 2,
    figure;
    if strcmp(model.type, 'cgpdisim'),
      numRows = 3;
      numCols = numGenes - 1;
      inputCol = floor((numCols+1)/2);
  
      subplot(numRows, numCols, inputCol + numCols);
    else
      numRows = 2;
      numCols = numGenes;
      inputCol = floor((numCols+1)/2);
      
      subplot(numRows, numCols, inputCol);
    end
  else
    figure(j);
    subplot(numPlots+1, totgenes, totgenes + thisgene);
  end
  bh = fill([predt', fliplr(predt')], [predF' + 2*sqrt(varF'), fliplr(predF' - 2*sqrt(varF'))], [.8 .8 .8], 'EdgeColor', [.8, .8, .8]);
  hold on,
  lin = plot(predt, predF, '-');
  hold off
  title(sprintf('Inferred %s protein', tf), 'fontsize', FONTSIZE);
  set(bh, 'lineWidth', LINEWIDTH);
  set(lin, 'lineWidth', 1.5*LINEWIDTH);
  %set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', FONTSIZE, 'xlim', [min(model.comp{j}.t) ...
                                                    max(model.comp{j} ...
                                                    .t)]);
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', []);
    axis tight
    v = axis;
    axis([v(1)-.05*(v(2)-v(1)), v(2)+.05*(v(2)-v(1)), ...
	  v(3)-.05*(v(4)-v(3)), v(4)+.05*(v(4)-v(3))])
%  set(gca, 'ylim', [-0.1 0.4]);
  
  for index = 1:numPlots
    if nargin >= 2,
      if index==1,
	subplot(numPlots+1, totgenes, thisgene);
      else
	subplot(numPlots+1, totgenes, index*totgenes+thisgene);
      end
    else
      if strcmp(model.type, 'cgpdisim'),
	if index == 1,
	  subplot(numRows, numCols, inputCol);
	else
	  subplot(numRows, numCols, (numRows-1)*numCols + index - 1);
	end
      else
	subplot(numRows, numCols, (numRows-1)*numCols + index);
      end
    end
    bh = fill([predt', fliplr(predt')], [predExprs(:,selectgenes(index))' + 2*sqrt(varExprs(:,selectgenes(index))'), fliplr(predExprs(:,selectgenes(index))' - 2*sqrt(varExprs(:,selectgenes(index))'))], [.8, .8, .8], 'EdgeColor', [.8, .8, .8]);
    hold on,
    lin = plot(predt, predExprs(:,selectgenes(index)), '-');
    lin = [lin plot(times, y{j}(:,selectgenes(index)), 'rx')];
    lin1 = errorbar(times, y{j}(:,selectgenes(index)), 2*sqrt(yvar{j}(:,selectgenes(index))), ...
                    'rx');
    hold off
    if strcmp(model.type, 'cgpdisim') && selectgenes(index) == 1,
      titleText = sprintf('%s mRNA (input)', tf);
%      lin = [lin plot(predt, predDi, 'm-')];     
    else
      titleText = sprintf('%s mRNA', genenames{selectgenes(index)});
    end
    title(titleText, 'fontsize', FONTSIZE);
    set(bh, 'lineWidth', LINEWIDTH);
    set(lin, 'lineWidth', 1.5*LINEWIDTH);
    set(lin, 'markersize', MARKERSIZE);
    set(lin1, 'lineWidth', LINEWIDTH);
    set(gca, 'fontname', 'arial', 'fontsize', FONTSIZE, 'xlim', [min(predt) ...
                        max(predt)]);
    set(gca, 'YTick', []);
    if nargin < 2,
      if strcmp(model.type, 'cgpdisim') && index == 1,
	set(gca, 'XTickLabel', []);
      else
	xlabel('Time (h)')
      end
    else
      if index ~= numPlots,
	set(gca, 'XTickLabel', []);
      else
	xlabel('Time (h)')
      end
    end
    axis tight
    v = axis;
    axis([v(1)-.05*(v(2)-v(1)), v(2)+.05*(v(2)-v(1)), ...
	  v(3)-.05*(v(4)-v(3)), v(4)+.05*(v(4)-v(3))])
  end
end
