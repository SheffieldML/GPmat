function gpsimMapEcoliResults(model, type, expNo, scale, printResults)

% GPSIMMAPECOLIRESULTS Plot the results from the MAP script.
% FORMAT
% DESC plots results from the MAP script, printing them if the
% given flag is set to true.
% ARG model : the model for which the results are being plotted.
% ARG type : gives a descriptive label for the type of the
% experiment.
% ARG expNo : number identifier for the experiment.
% ARG scale : the scaling factor applied to normalise. 
% ARG printResults : function prints the plots if set to true
% (default is false).
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008
%
% SEEALSO : demEcoliMap1

% SHEFFIELDML

type(1) = upper(type(1));
switch nargin
 case 3
  scale = 1;
  printResults = false;
 case 4
  printResults = false;
end

gene = {'dinF', 'dinI', 'lexA', 'recA', 'recN', 'ruvA', 'ruvB', ...
                  'sbmC', 'sulA', 'umuC', 'umuD', 'uvrB', 'yebG', 'yjiW' };
if model.comp{1}.numGenes == 8
  targetInd = [2 3 6 7 9 11 12 14];
  gene = gene(targetInd);
end  
  
predt = model.comp{1}.mapt;
predExprsFull = [];
for i = 1:length(model.comp)
    predExprsFull = [predExprsFull; model.comp{i}.ypred];
end
meanPredExprs = ones(size(predExprsFull, 1), 1)*mean(predExprsFull);
scalePredExprs = sqrt(var(predExprsFull - meanPredExprs));
scaleMat = ones(size(predExprsFull, 1), 1)*scalePredExprs;
% predExprsFull = meanPredExprs + (predExprsFull-meanPredExprs)./scaleMat;
startVal = 1;
endVal = length(predt);
for i = 1:length(model.comp)
    predExprs{i} = predExprsFull(startVal:endVal, :);
    startVal = startVal + length(predt);
    endVal = endVal + length(predt);
end

modelB = model.comp{1}.B.*scale;
modelS = model.comp{1}.S.*scale;
modelAlpha = model.comp{1}.alpha.*scale;

for j = 1:length(model.comp)
  
  scalePred = sqrt(var(model.comp{j}.f));
  truef = [0.16 0 0.19 0.75 1];
%  truef = truef/sqrt(var(truef))*scalePred;
    
  % Figure 2(a) histograms;
%   B = trueParams.B./scale; % From Martino paper ... but don't know the
%                            % scale
%   B = B/mean(B)*mean(modelB); % do a rough rescaling so
%                                         % that the scales match.
%   S = trueParams.S./scale; 
%   S = S/mean(S)*mean(modelS); 
%   D = trueParams.D; % From Martino paper, again
%                                       % we know the scale because
%                                       % p21 is fixed to 0.8.  
  
  % Martino f from Figure 2(b), again measured with a ruler.
  if model.comp{j}.ngParam
    nonLinearity = 'multi';
  else
    nonLinearity = model.comp{j}.nonLinearity;
  end
  switch nonLinearity
   case 'linear'
    figure, lin = plot(model.comp{j}.mapt,model.comp{j}.f, '-');
    hold on;
    bh = plot(model.comp{j}.mapt, model.comp{j}.f+2*sqrt(model.comp{j}.varF),'--');
    bh =[bh plot(model.comp{j}.mapt, model.comp{j}.f-2*sqrt(model.comp{j}.varF),'--')];
    lin = [lin plot(model.comp{j}.t(1:5), truef, 'rx')];
   case 'multi'
    figure, 
    predF = exp(model.comp{j}.f);
    varFplus = exp(model.comp{j}.f+2*sqrt(model.comp{j}.varf));
    varFminus = exp(model.comp{j}.f-2*sqrt(model.comp{j}.varf));
    lin = plot(model.comp{j}.mapt, predF, '-');
    hold on;
    bh = plot(model.comp{j}.mapt, varFplus, '--');
    bh = [bh plot(model.comp{j}.mapt, varFminus, '--')];
%    truef = truef/sqrt(var(truef))*sqrt(var(predF));
    title('Inferred LexA Activity', 'fontsize', 20);
%    lin = [lin plot(model.comp{j}.t(1:5), truef, 'rx')];    
   otherwise
     figure,
% 
%     predF = exp(model.comp{j}.f);
%     varFplus = exp(model.comp{j}.f+2*sqrt(model.comp{j}.varf));
%     varFminus = exp(model.comp{j}.f-2*sqrt(model.comp{j}.varf));
%     lin = plot(model.comp{j}.mapt, predF-min(predF), '-');
%     hold on;
%     bh = plot(model.comp{j}.mapt, varFplus-min(predF), '--');
%     bh = [bh plot(model.comp{j}.mapt, varFminus-min(predF), '--')];
%     truef = truef/sqrt(var(truef))*sqrt(var(predF));   
%     
    
    subplot(2,1,1);
    lin = plot(model.comp{j}.mapt, model.comp{j}.f, '-');
    hold on;
    bh = plot(model.comp{j}.mapt, model.comp{j}.f+2*sqrt(model.comp{j}.varf),'--');
    bh = [bh plot(model.comp{j}.mapt, model.comp{j}.f-2* ...
                                           sqrt(model.comp{j}.varf),'--')];
    lin = [lin plot(model.comp{j}.t(1:5), truef, 'rx')];
    title('f(t) Profile');
    
    func = str2func(model.comp{j}.nonLinearity);
    subplot(2,1,2); lin = [lin plot(model.comp{j}.mapt, func(model.comp{j}.f), '-')];
    hold on;
    bh = [bh plot(model.comp{j}.mapt, func(model.comp{j}.f+2*sqrt(model.comp{j}.varf)),'--')];
    bh = [bh plot(model.comp{j}.mapt, func(model.comp{j}.f-2* ...
                                           sqrt(model.comp{j}.varf)),'--')];
    lin = [lin plot(model.comp{j}.t(1:5), func(truef), 'rx')];
    title('g(f) Profile');
  end

  set(bh, 'lineWidth', 2);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
 
  kernName = model.comp{j}.kern.type;
  kernName(1) = upper(kernName(1));
  
  switch nonLinearity
   case 'linear'
    set(gca, 'ylim', [-2 4])
   case {'exp', 'negLogLogit'}
    set(gca, 'ylim', [0 8])
%    otherwise
%     set(gca, 'ylim', [0 2])
  end
 
  fileName = ['dem' type num2str(expNo) kernName ...
              nonLinearity '_profile' num2str(j) '_slide'];
  if printResults 
    print('-deps', ['./results/' fileName]);
  end
  pos = get(gcf, 'paperposition');
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  if printResults
    print('-dpng', ['./results/' fileName])
  end
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos)
  switch nonLinearity
   case 'linear'
    xStart = min(model.comp{1}.mapt);
    xEnd = max(model.comp{1}.mapt);  
    axis([xStart xEnd -2 3.5]);
  end
  
  
  
  for index = 1:model.comp{j}.numGenes
      figure;
      lin = plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index)*scale(index), '-');
      hold on;
      bh = plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index)*scale(index) + ...
                2*sqrt(model.comp{j}.ypredVar(:,index)*scale(index)*scale(index)),'--');
      bh = [bh plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index)*scale(index) - ...
                2*sqrt(model.comp{j}.ypredVar(:,index)*scale(index)*scale(index)) ,'--')];      
      lin = [lin plot(model.comp{j}.t, model.comp{j}.y(:,index)*scale(index), 'rx')];
      titleText = [gene{index} ' mRNA'];
      title(titleText,'fontsize', 20);
      
      textB = ['B = ' num2str(modelB(index))];
      textD = ['D = ' num2str(model.comp{j}.D(index))];
      textS = ['S = ' num2str(modelS(index))];
      textAlpha = ['Alpha = ' num2str(modelAlpha(index))];
      textGamma = ['Gamma = ' num2str(model.comp{j}.gParam(1,index))];
      texPos = ylim;
      scalePos = (texPos(2) - texPos(1))/14;
      texh = text(35, texPos(1)+5*scalePos, textB);
      texh = [texh text(35, texPos(1)+4*scalePos,textD)];
      texh = [texh text(35, texPos(1)+3*scalePos, textS)];  
      texh = [texh text(35, texPos(1)+2*scalePos, textAlpha)];        
      texh = [texh text(35, texPos(1)+1*scalePos, textGamma)];        
      set(bh, 'lineWidth', 3);
      set(lin, 'lineWidth', 4);
      set(lin, 'markersize', 20);
      set(texh, 'fontsize', 16);
      set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [min(model.comp{j}.mapt) max(model.comp{j}.mapt)])

      if printResults
          fileName = ['dem' type num2str(expNo) kernName ...
              nonLinearity '_ExprsProfile_Rep' num2str(j) '_Gene' num2str(index)];
          print('-deps', ['./results/' fileName]);
          pos = get(gcf, 'paperposition');
          origpos = pos;
          pos(3) = pos(3);
          pos(4) = pos(4);
          set(gcf, 'paperposition', pos);
          lineWidth = get(gca, 'lineWidth');
          set(gca, 'lineWidth', lineWidth*2);
          print('-dpng', ['./results/' fileName]);
          set(gca, 'lineWidth', lineWidth);
          set(gcf, 'paperposition', origpos);
      end
  end  
end

order = 1:model.comp{1}.numGenes;

counter = 0;

% Plot first basal transcription rates.
figure
bar([modelB(order)]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', gene)
if printResults
  fileName = ['dem' type num2str(expNo) kernName ...
              nonLinearity '_basal'];
  print('-deps', ['./results/' fileName]);
end
pos = get(gcf, 'paperposition');
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
if printResults
  print('-dpng', ['./results/' fileName]);
end
set(gcf, 'paperposition', origpos);
set(gca, 'lineWidth', lineWidth);

% Plot the sensitivities.
figure
bar([modelS(order)]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', gene);
if printResults
  fileName = ['dem' type num2str(expNo) kernName ...
              nonLinearity '_sensitivity'];
  print('-deps', ['./results/' fileName]);
end
pos = get(gcf, 'paperposition');
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
if printResults
  print('-dpng', ['./results/' fileName]);
end
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);

% Finally plot degradation rates.
figure
bar([model.comp{1}.D(order)]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', gene);
if printResults
  fileName = ['dem' type num2str(expNo) kernName ...
              nonLinearity '_decay'];
  print('-deps', ['./results/' fileName]);
end
pos = get(gcf, 'paperposition');
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
if printResults
  print('-dpng', ['./results/' fileName]);
end
set(gcf, 'paperposition', origpos);
set(gca, 'lineWidth', lineWidth);


