% DEMTOYPROBLEM2 Display results from toy problem point at a time.

% SHEFFIELDML

load demToyProblem1.mat

bw = true; %false;

predt = [linspace(0, 18, 100)]';
figure(1), clf
order = [3 4 6 2 5];
lin1 = plot(t, truey(:, order-1));
hold on
lin2 = [];
counter = 0;
presentOrder = 0;
for i = order
  presentOrder = presentOrder + 1;
  for j = 1:length(model.timesCell{i})
    counter = counter + 1;
    indices(counter, 1) = i;
    indices(counter, 2) = j;
    indices(counter, 3) = presentOrder;
  end
end
ind = [1];
tvals = cell(length(model.timesCell));
for i = 1:size(indices, 1)
  offset = 0;
  for j = 1:indices(i, 1)-1
    offset = offset + length(model.timesCell{j});
  end
  ind = [ind indices(i, 2)+offset];
  tvals{indices(i, 1)} = [tvals{indices(i, 1)}; model.timesCell{indices(i, 1)}(indices(i, ...
                                                    2))];
  
  proteinKern = model.kern.comp{1};
  K = rbfKernCompute(proteinKern, 0, predt);
  counter = 0;
  for j=order
    counter = counter + 1;
    if ~isempty(tvals{j})
      K = [K; real(simXrbfKernCompute(model.kern.comp{j}, proteinKern, ...
                                      tvals{j}, predt))];
    end
  end  
  invK = pdinv(model.K(ind, ind));
  obsY = model.m(ind, 1);
  predF = K'*invK*obsY;
  varF = kernDiagCompute(proteinKern, predt) - sum(K.*(invK*K), 1)';
  
  figure(1)
  lin2 = [ plot(model.timesCell{indices(i, 1)}(indices(i, 2)), ...
               [repmat(NaN, 1, indices(i, 3)-1) model.y(ind(end) - 1)], '.')];
  set(lin1, 'lineWidth', 2);
  set(lin1, 'markersize', 20);
  set(lin2, 'lineWidth', 4);
  set(lin2, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  if bw
    fileName = ['demToyProblem2_genes'];
    set(lin2, 'color', [0 0 0]);
    print('-deps', ['../tex/diagrams/' fileName num2str(i)]);
  else
    fileName = ['demToyProblem2bw_genes'];
    print('-depsc', ['../tex/diagrams/' fileName num2str(i)]);
    pos = get(gcf, 'paperposition');
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
  
  

  
  
  figure(2), clf  
  hold on
  stdVals = sqrt(varF);
  fillColor = [0.7 0.7 0.7];
  fill([predt; predt(end:-1:1)], ...
       [predF; predF(end:-1:1)] ...
       + 2*[stdVals; -stdVals(end:-1:1)], ...
       fillColor,'EdgeColor',fillColor)
  lin = plot(t, truef, 'r-');  
  lin = [lin plot(predt, predF, '-')];
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  set(gca, 'ylim', [-2 4])
  if bw
    fileName = ['demToyProblem2bw_infered'];
    print('-deps', ['../tex/diagrams/' fileName num2str(i)]);
    set(lin, 'color', [0 0 0]);
    set(lin, 'color', [0 0 0]);
  else
    fileName = ['demToyProblem2_infered'];
    print('-depsc', ['../tex/diagrams/' fileName num2str(i)]);
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    print('-dpng', ['../html/' fileName])
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos)
  end
end


%save demToyProblem2.mat
