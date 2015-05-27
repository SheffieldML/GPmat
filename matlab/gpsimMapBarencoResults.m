function gpsimMapBarencoResults(model, type, expNo, printResults, scale)

% GPSIMMAPBARENCORESULTS Plot the results from the MAP script.
% FORMAT
% DESC plots results from the MAP script, printing them if the
% given flag is set to true.
% ARG model : the model for which the results are being plotted.
% ARG type : gives a descriptive label for the type of the
% experiment.
% ARG expNo : number identifier for the experiment.
% ARG printResults : function prints the plots if set to true
% (default is false).
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006
%
% MODIFIED : Pei Gao, 2008
%
% SEEALSO : demBarencoMap1

% SHEFFIELDML

for j = 1:length(model.comp)
  model.comp{j} = gpsimMapUpdateYpredVar(model.comp{j});
end
  
type(1) = upper(type(1));
switch nargin
 case 3
  scale = 1;
  printResults = false;
 case 4
  scale = 1;
end

modelB = model.comp{1}.B.*scale;
if model.comp{1}.ngParam
  modelS = (model.comp{1}.S./(1+model.comp{1}.gParam)).*scale;
else
  modelS = model.comp{1}.S.*scale;
end
modelS = modelS/modelS(4);

if model.comp{1}.ngParam
  nonLinearity = 'multi';
else
  nonLinearity = model.comp{1}.nonLinearity;
end

for j = 1:length(model.comp)
  
  if model.comp{1}.ngParam
    scalePred = sqrt(var(exp(model.comp{j}.f)));
  else
    scalePred = sqrt(var(model.comp{j}.f)); 
  end
  
  % Info from Martino Paper:
  % 'True f' from Figure 3.
  % Got these figures with a ruler ...
  % Don't actually plot these below, but they are stored for reference
  truef = [0 1.6 2.6 2.5 2.6 1.6 0.9];
  truef = truef/sqrt(var(truef))*scalePred;
  
  
  % Figure 2(a) histograms;
  B = [2.6 1.5 0.5 0.2 1.35]; % From Martino paper ... but don't know the scale
  B = B/mean(B)*mean(modelB); % do a rough rescaling so
                                       % that the scales match.
  S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                               % know the scale, because p21 is
                               % fixed to 1.
  D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.  
  
  % Martino f from Figure 2(b), again measured with a ruler.
  barencof = [0.0000000 200.5201100 355.5216125 205.7574913 135.0911372 ...
              145.1080997 130.7046969; 0.0000000 184.0994134 308.4759200 ...
              232.1775328 153.6595161 85.7272235  168.0910562; 0.0000000 ...
              230.2262511 337.5994811 276.9416540 164.5044287 127.8653452 173.6112139];
  barencof = barencof/(1.8*mean(S))*mean(modelS);
  barencof = barencof./(sqrt(var(barencof,0,2))*ones(1,7))*scalePred;
  
  switch nonLinearity
   case 'linear'
    figure, lin = plot(model.comp{j}.mapt,model.comp{j}.f, '-');
    hold on;
    bh = plot(model.comp{j}.mapt, model.comp{j}.f+2*sqrt(model.comp{j}.varf),'--');
    bh =[bh plot(model.comp{j}.mapt, model.comp{j}.f-2* ...
                 sqrt(model.comp{j}.varf),'--')];
   case 'multi'
    figure, 
%     subplot(1,2,1); 
%     lin = plot(model.comp{j}.mapt,model.comp{j}.f, '-');
%     hold on;
%     bh = plot(model.comp{j}.mapt, model.comp{j}.f+2*sqrt(model.comp{j}.varf),'--');
%     bh =[bh plot(model.comp{j}.mapt, model.comp{j}.f-2* ...
%                  sqrt(model.comp{j}.varf),'--')];    
%     
%     subplot(1,2,2);
    lin = plot(model.comp{j}.mapt, exp(model.comp{j}.f), '-');
    hold on;
    bh = plot(model.comp{j}.mapt, exp(model.comp{j}.f+2*sqrt(model.comp{j}.varf)),'--');
    bh = [bh plot(model.comp{j}.mapt, exp(model.comp{j}.f-2* ...
                                          sqrt(model.comp{j}.varf)),'--')];
   otherwise
    func = str2func(model.comp{j}.nonLinearity);
    figure, lin = plot(model.comp{j}.mapt, func(model.comp{j}.f), '-');
    hold on;
    bh = plot(model.comp{j}.mapt, func(model.comp{j}.f+2*sqrt(model.comp{j}.varf)),'--');
    bh = [bh plot(model.comp{j}.mapt, func(model.comp{j}.f-2*sqrt(model.comp{j}.varf)),'--')];
  end
%  lin = [lin plot(0:2:12, barencof(j,:),'rx')];
  title('Inferred p53 protein','fontsize', 20);
  set(bh, 'lineWidth', 3);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [0 12])
  kernName = model.comp{j}.kern.type;
  kernName(1) = upper(kernName(1));
  
  switch nonLinearity
   case 'linear'
    set(gca, 'ylim', [-2 4])
   case {'exp', 'negLogLogit'}
    set(gca, 'ylim', [0 4])
  end

  fileName = ['demBarenco' type num2str(expNo) kernName ...
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
    lin = plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index), '-');
    hold on;
    bh = plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index)+2* ...
              sqrt(model.comp{j}.ypredVar(:,index)),'--');
    bh = [bh plot(model.comp{j}.mapt, model.comp{j}.ypred(:,index)-2* ...
                  sqrt(model.comp{j}.ypredVar(:,index)),'--')];
    lin = [lin plot(model.comp{j}.t, model.comp{j}.y(:,index), 'rx')];
    lin1 = errorbar(model.comp{j}.t, model.comp{j}.y(:,index), 2* ...
                    sqrt(model.comp{j}.yvar(:,index)), 'rx');
    titleText = [num2str(index) '-th Gene'];
    title(titleText);
    set(bh, 'lineWidth', 3);
    set(lin, 'lineWidth', 4);
    set(lin, 'markersize', 20);
    set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', ...
             [min(model.comp{j}.mapt) max(model.comp{j}.mapt)])
    
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


order = [1 5 3 4 2];

counter = 0;

% Plot first basal transcription rates.
figure
bar([modelB(order); B]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', 'BIK'})
if printResults
  fileName = ['demBarenco' type num2str(expNo) kernName ...
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
  print('-dpng', ['./results/' fileName])
end
set(gcf, 'paperposition', origpos);
set(gca, 'lineWidth', lineWidth);

% Plot the sensitivities.
figure
bar([modelS(order); S]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
if printResults
  fileName = ['demBarenco' type num2str(expNo) kernName ...
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
  print('-dpng', ['./results/' fileName])
end
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);

% Finally plot degradation rates.
figure
bar([model.comp{1}.D(order); D]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
if printResults
  fileName = ['demBarenco' type num2str(expNo) kernName ...
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
  print('-dpng', ['./results/' fileName])
end
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);
