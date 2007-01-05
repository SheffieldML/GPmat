function gpsimMapBarencoResults(model, type, expNo, printResults)

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
% SEEALSO : demBarencoMap1

% GPSIM

type(1) = upper(type(1));
if nargin < 4
  printResults = false;
end
for j = 1:length(model.comp)
  
  %/~
  %   % Generate predictions of the functions.
  %   % to do this we need to compute the K_xf portions of the kernel
  %   % (simXrbfKernCompute does this for us).
  %   predt = [linspace(-2, 14, 161) 0:2:12]';
  %   proteinKern = kernCreate(origModel.comp{1}.t, 'rbf');
%   proteinKern.inverseWidth = ...
%       origModel.comp{j}.kern.comp{1}.inverseWidth;
%   K = [];
%   for i=1:origModel.comp{j}.kern.numBlocks
%     K = [K; simXrbfKernCompute(origModel.comp{j}.kern.comp{i}, proteinKern, ...
%                                origModel.comp{j}.t, predt)];
    
%   end
%   predF = K'*origModel.comp{j}.invK*origModel.comp{j}.y;
%   varF = kernDiagCompute(proteinKern, predt) - sum(K.*(origModel.comp{j}.invK*K), 1)';
  
%   % Take out predictions at data points.
%   % Use them to get the scale for the other data.
%   dataF = predF(end-6:end);
%   dataVarF = varF(end-6:end);
%   predF(end-6:end) = [];
%   varF(end-6:end) = [];
%   predt(end-6:end) = [];
%   scalePred = sqrt(var(dataF));
%~/  

  scalePred = sqrt(var(model.comp{j}.g));

  % Info from Martino Paper:
  % 'True f' from Figure 3.
  % Got these figures with a ruler ...
  % Don't actually plot these below, but they are stored for reference
  truef = [0 1.6 2.6 2.5 2.6 1.6 0.9];
  truef = truef/sqrt(var(truef))*scalePred;
  
  
  % Figure 2(a) histograms;
%  B = [2.6 1.5 0.5 0.2 1.35]; % From Martino paper ... but don't know the scale
%  B = B/mean(B)*mean(model.comp{j}.B); % do a rough rescaling so
                                       % that the scales match.
  S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                               % know the scale, because p21 is
                               % fixed to 1.
  D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.  
  
  % Martino f from Figure 2(b), again measured with a ruler.
  barencof = [0 2.7 3.9 2.3 1.5 1.6 1.4]/(1.8*mean(S))*mean(model.comp{j}.S);
  barencof = barencof/sqrt(var(barencof))*scalePred;
  
  
  switch model.comp{j}.nonLinearity
   case 'linear'
    figure, lin = plot(model.comp{j}.mapt,model.comp{j}.f, '-');
    hold on;
    bh = plot(model.comp{j}.mapt, model.comp{j}.f+2*sqrt(model.comp{j}.varf),'--');
    bh =[bh plot(model.comp{j}.mapt, model.comp{j}.f-2*sqrt(model.comp{j}.varf),'--')];
   otherwise
    func = str2func(model.comp{j}.nonLinearity);
    figure, lin = plot(model.comp{j}.mapt, func(model.comp{j}.f), '-');
    hold on;
    bh = plot(model.comp{j}.mapt, func(model.comp{j}.f+2*sqrt(model.comp{j}.varf)),'--');
    bh = [bh plot(model.comp{j}.mapt, func(model.comp{j}.f-2*sqrt(model.comp{j}.varf)),'--')];
  end
  lin = [lin plot(0:2:12, barencof,'rx')];
  set(bh, 'lineWidth', 3);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  kernName = model.comp{j}.kern.type;
  kernName(1) = upper(kernName(1));
  
  switch model.comp{j}.nonLinearity
   case 'linear'
    set(gca, 'ylim', [-2 4])
   case {'exp', 'negLogLogit'}
    set(gca, 'ylim', [0 6])
  end
  nonLinearity = model.comp{j}.nonLinearity;
  nonLinearity(1) = upper(nonLinearity(1));
  fileName = ['demBarenco' type num2str(expNo) kernName ...
              nonLinearity '_profile' num2str(j) '_slide'];
  if printResults 
    print('-depsc', ['../tex/diagrams/' fileName]);
  end
  pos = get(gcf, 'paperposition');
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  if printResults
    print('-dpng', ['../html/' fileName])
  end
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos)
  switch model.comp{j}.nonLinearity
   case 'linear'
    axis([-2 14 -0.5 5.5]);
  end
end