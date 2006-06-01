% DEMBARENCOMAP1 Optimise model using MAP approximation.

% GPSIM

% Map learning of f given replicated data and known parameters
% MAP-Laplace optimisation of the length-scale
% Constraining TF concentration to be positive

% Run demBarenco1.m first, to get data and initial B,S,D and lengthscale

%data collection times
times = [0 2 4 6 8 10 12]';

%spacing used in part of Neil's code pasted at the bottom of this file
predt = linspace(-2, 14, 161);

%spacing used for integration here (I've made it the same, but can be finer for greater accuracy)

load demBarenco1.mat

left=-2;
right=14;
step = 0.1;
mapt=[];
mapt=[mapt,left:step:times(1)];
times_index = [length(mapt)];
for i=1:length(times)-1
  mapt=[mapt,times(i)+step:step:times(i+1)];
  times_index = [times_index,length(mapt)];
end;
mapt=[mapt,times(length(times))+step:step:right]'; %time vector
npts = length(mapt);

%gam=1; %Initialise gamma - inverse length scale.

%model.K = zeros(npts,npts);
%for i=1:npts
%   for j=1:npts
%        model.K(i,j) = exp(-gam*(mapt(i)-mapt(j))^2/2);
%    end;
%end;

%These come from Neil's results
origModel = model;
ngenes = origModel.comp{1}.numGenes;
B = origModel.comp{1}.B;
D = origModel.comp{1}.D;
S = origModel.comp{1}.S;


nreps = length(y); %number of replicates

fstart=zeros(nreps,length(mapt)); %keep initial state for viewing progress

%initialisations
f = fstart;                     %TF concentrations for each replicate
varf=zeros(nreps,length(mapt));   %Variances
W = zeros(length(f),length(f)); %Hessian of -ve log-likelihood
df = zeros(1,length(mapt));            %gradient of log-likelihood
dg = cell(1,nreps);

options = gpOptions('ftc');
options.kern = kernCreate(mapt, 'mlp');
model = gpCreate(1, 1, mapt, f(1, :)', options);

model.kern.weightVariance = 30;
model.kern.biasVariance = 1000;

x = zeros(ngenes,npts); %numerical estimate of x(t)

gamvec{i}=kernExtractParam(model.kern); %vector of gamma estimates

prev_err = ones(nreps,1);

iters = 70; %Number of optimisation iterations
for ii=1:iters %Start optimisation
    
  %  %Update the kernel   
  model.K = kernCompute(model.kern, mapt);
  model.Kinv = pdinv(model.K);
  %  model.K=zeros(npts,npts);
  %  for i=1:npts
  %    for j=1:npts
  %      model.K(i,j) = exp(-gam*(mapt(i)-mapt(j))^2/2);
  %    end;
  %  end;   
  
  for rep=1:nreps %Work out likelihood gradient for each replicate
    
    x_data = y{rep};    %data
    x_var = yvar{rep};  %variance
    
    %Run the dynamics forward from t=0
    for i=1:ngenes
      x(i,times_index(1)) = B(i)/D(i);
      integral=0;
      for j=times_index(1)+1:npts    
        tfs=mapt(j);
        integral=integral+exp(f(rep,j)+D(i)*tfs)*step;
        x(i,j) = B(i)/D(i) + S(i)*exp(-D(i)*tfs)*integral;
      end;
    end;
    
    % Functional gradient of log-likelihood
    for j=times_index(1):npts 
      temp=0;
      for k=1:length(times)
        arg = times(k)-mapt(j);
        if arg > 0
          for i=1:ngenes
            lambda=1/x_var(k,i);
            factor=(x(i,times_index(k))-x_data(k,i))*lambda;
            temp=temp+factor*exp(f(rep,j)-D(i)*arg)*S(i);
          end;
        end;
      end;
      df(j) = -temp*step;
    end;
    
    % Hessian of -ve log-likelihood
    for j=times_index(1):npts 
      for l=times_index(1):npts 
        temp=0;
        for k=1:length(times)
          arg1 = times(k)-mapt(j);
          arg2 = times(k)-mapt(l);
          if arg1 > 0 && arg2 > 0
            for i=1:ngenes
              lambda=1/x_var(k,i);
              temp=temp+lambda*exp(f(rep,j)+f(rep,l)-D(i)*(arg1+arg2))*S(i)^2;
            end;
          end;
        end;
        W(j,l) = temp;
      end;
    end;
    W=W*step^2;
    prefac = pdinv(model.Kinv+W);
    
    fvec = f(rep,:);
    
    %if abs(prev_err(rep))<1
    %    f(rep,:) = (prefac*(W*fvec' + df'))';   %Full Newton iteration
    %else
    %    
    
    new_drn = (prefac*(df' - model.Kinv*fvec'))';
    f(rep,:) = f(rep,:) + 0.5*new_drn;  %Partial Newton iteration
    
    %end;
    varf(rep,:) = diag(prefac)';
    
    prev_err(rep)=sum(fvec - f(rep,:))
    fvec = f(rep,:);
    
    if sum(abs(prev_err)<0.01)
      updateKernParam=true;
    else
      updateKernParam=false;
    end
    
    %Find gradient with respect to kernel parameters if updating gamma
    
    if updateKernParam
      %dKdg=zeros(npts,npts);
      %for i=1:npts
      %  for j=1:npts
      %    dKdg(i,j) = -0.5*(mapt(i)-mapt(j))^2*model.K(i,j);
      %  end;
      %end;   
      %This part is same as linear case
      %dg{rep} = 0.5*fvec*model.Kinv*dKdg*model.Kinv*fvec'-0.5*trace(pdinv(pdinv(W)+model.K)*dKdg);
      %This is the additional part due to W's f dependence
      %term1 = -diag(prefac*W);
      %term2 = prefac*model.Kinv*dKdg*df';
      %dg{rep} = dg{rep} + term1'*term2;
      dg{rep} = kernGradient(model.kern, mapt, ...
                         0.5*model.Kinv*fvec'*fvec*model.Kinv ...
                         - 0.5*pdinv(pdinv(W)+model.K) ...
                         -model.Kinv*prefac*diag(prefac*W)*df);
    end 
  end; %replicates
  
  %Update kernel parameters by simple gradient ascent
  
  if updateKernParam
    eta=0.02;
    param = kernExtractParam(model.kern);
    for i = 1:length(dg)
      param(1:end-1) = param(1:end-1) + eta*sum(dg{i}(1:end-1));
    end
    model.kern = kernExpandParam(model.kern, param);
    gamvec{end+1}=kernExtractParam(model.kern);
  end 
  
  %Plot convergence for 1st replicate
  if mod(ii,1)==0
    figure(1);
    expMean = exp(f + varf/2);  % Logs are normally distributed
                                     % ... recover mean in exp space.
    expVar = (exp(varf)-1).*exp(2*f + varf); % Logs are
                                                  % normally
                                                  % distributed
                                                  % ... recover
                                                  % variance in exp
                                                  % space.

    plot(mapt, expMean(1, :));
    hold on
    plot(mapt, expMean(1, :)+2*sqrt(expVar(1, :)));
    plot(mapt, expMean(1, :)-2*sqrt(expVar(1, :)));
%    plot(mapt,exp(f(1,:)));
%    hold on;
%    plot(mapt,exp(f(1,:)+2*sqrt(varf(1,:))));
%    plot(mapt,exp(f(1,:)-2*sqrt(varf(1,:))));
    hold off;
    pause(0.1);
  end % mod
  
end %iters

%Plot results 

for j = 1:3
  
  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).
  predt = [linspace(-2, 14, 161) 0:2:12]';
  proteinKern = kernCreate(origModel.comp{1}.t, 'rbf');
  proteinKern.inverseWidth = ...
      origModel.comp{j}.kern.comp{1}.inverseWidth;
  K = [];
  for i=1:origModel.comp{j}.kern.numBlocks
    K = [K; simXrbfKernCompute(origModel.comp{j}.kern.comp{i}, proteinKern, ...
                               origModel.comp{j}.t, predt)];
    
  end
  predF = K'*origModel.comp{j}.invK*origModel.comp{j}.y;
  varF = kernDiagCompute(proteinKern, predt) - sum(K.*(origModel.comp{j}.invK*K), 1)';
  
  % Take out predictions at data points.
  % Use them to get the scale for the other data.
  dataF = predF(end-6:end);
  dataVarF = varF(end-6:end);
  predF(end-6:end) = [];
  varF(end-6:end) = [];
  predt(end-6:end) = [];
  scalePred = sqrt(var(dataF));
  
  % Info from Martino Paper:
  % 'True f' from Figure 3.
  % Got these figures with a ruler ...
  % Don't actually plot these below, but they are stored for reference
  truef = [0 1.6 2.6 2.5 2.6 1.6 0.9];
  truef = truef/sqrt(var(truef))*scalePred;
  
  
  % Figure 2(a) histograms;
  B = [2.6 1.5 0.5 0.2 1.35]; % From Martino paper ... but don't know the scale
  B = B/mean(B)*mean(origModel.comp{1}.B); % do a rough rescaling so
                                       % that the scales match.
  S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                               % know the scale, because p21 is
                               % fixed to 1.
  D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.
  
  % Martino f from Figure 2(b), again measured with a ruler.
  barencof = [0 2.7 3.9 2.3 1.5 1.6 1.4]/(1.8*mean(S))*mean(origModel.comp{1}.S);
  barencof = barencof/sqrt(var(barencof))*scalePred;
  
  
  figure, lin = plot(mapt, expMean(j, :), '-');
  hold on
  plot(mapt, expMean(j, :)+2*sqrt(expVar(j, :)), '--');
  plot(mapt, expMean(j, :)-2*sqrt(expVar(j, :)), '--');
%  figure, lin = plot(mapt,exp(f(j,:)), '-');
%  hold on;
%  plot(mapt, exp(f(j,:)+2*sqrt(varf(j,:))),'--')
%  plot(mapt, exp(f(j,:)-2*sqrt(varf(j,:))),'--')
  lin = [lin plot(0:2:12, barencof,'x')];
  set(lin, 'lineWidth', 2);
  set(lin, 'markersize', 10);
  fileName = ['demBarenco1_profile' num2str(j)];
  %print('-deps', ['../tex/diagrams/' fileName]);
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
  axis([-2 14 -0.5 5.5]);
  
end




