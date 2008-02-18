% DEMTOYPROBLEM2 Display results from toy problem point at a time.

% GPSIM

load demToyProblem1.mat

predt = [linspace(0, 14, 100)]';

count = 0;
blockInd = [];
for i = 5;
  blockInd = [blockInd i];
  for j = 1:length(model.t);
    count = count + 1;
    if i<5
      yInd{count} = yInd{count-1};
    else
      if j == 1
        yInd{count} = j;
      else
        yInd{count} = [yInd{count-1} j];
      end
    end
    blockInc{count} = blockInd;
  end
end
for i = 4:-1:1
  blockInd = [blockInd i];
  for j = length(model.t);
    count = count + 1;
    if i<5
      yInd{count} = yInd{count-1};
    else
      if j == 1
        yInd{count} = j;
      else
        yInd{count} = [yInd{count-1} j];
      end
    end
    blockInc{count} = blockInd;
  end
end
  
for diaNo = 1:length(blockInc);
  ind = [];
  for i=blockInc{diaNo}
    ind = [ind yInd{diaNo}+(i-1)*length(model.t)];
  end
  tvals = model.t(yInd{diaNo}, 1);
  
  proteinKern = kernCreate(tvals, 'rbf');
  proteinKern.inverseWidth = ...
      model.kern.comp{1}.inverseWidth;
  K = [];
  for i=blockInc{diaNo}
    K = [K; simXrbfKernCompute(model.kern.comp{i}, proteinKern, ...
                               tvals, predt)];
    
  end
  invK = pdinv(model.K(ind, ind));
  obsY = model.y(ind, 1);
  startInd = 1;
  for i = 1:length(blockInc{diaNo})
    endInd = i*length(tvals);
    obsY(startInd:endInd) = obsY(startInd:endInd)-model.mu(blockInc{diaNo}(i));
    startInd = endInd + 1;
  end
  
  predF = K'*invK*obsY;
  varF = kernDiagCompute(proteinKern, predt) - sum(K.*(invK*K), 1)';
  
  % Take out predictions at data points.
  % Use them to get the scale for the other data.
  if diaNo == 1
%    dataF = predF(end-6:end);
%    dataVarF = varF(end-6:end);
    predF(end-6:end) = [];
    varF(end-6:end) = [];
    predt(end-6:end) = [];
    scalePred = sqrt(var(dataF));
  end

  figure(1), clf
  lin1 = plot(t, truey);
  hold on
  lin2 = [plot(tvals, reshape(model.y(ind, 1), length(tvals), length(blockInc{diaNo})), 'x')]
  set(lin1, 'lineWidth', 2);
  set(lin1, 'markersize', 20);
  set(lin2, 'lineWidth', 4);
  set(lin2, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  fileName = ['demToyProblem2_genes'];
  print('-depsc', ['../tex/diagrams/' fileName num2str(diaNo)]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  %print('-dpng', ['../html/' fileName])
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos)
  
  

  
  
  figure(2), clf, lin = plot(t, truef, 'r-');
  
  hold on,
  lin = [lin plot(predt, predF, '-')];
  bh = plot(predt, predF + 2*sqrt(varF), '--');
  bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
  set(bh, 'lineWidth', 3);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  fileName = ['demToyProblem2_infered'];
  set(gca, 'ylim', [-2 4])
  print('-depsc', ['../tex/diagrams/' fileName num2str(diaNo)]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  %print('-dpng', ['../html/' fileName])
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos)
end


%save demToyProblem2.mat