function model = gpsimMapFunctionalUpdateW(model)

% GPSIMMAPFUNCTIONALUPDATEW Update the data component of the Hessian.
% FORMAT
% DESC updates the portion of the Hessian associated with the
% likelihood of the data.
% ARG model : the model for which the data component of the Hessian is
% to be updated.
% RETURN model : the model with the data component of the Hessian up
% to date.
%
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalGradient, gpsimMapUpdatePosteriorCovariance
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

if isfield(model, 'includeNoise') && model.includeNoise
  noiseMat = ones(length(model.t),1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
else
  yvar = model.yvar;
end

if model.updateW
    intPoints = model.times_index(1)+1:(model.numMapPts);
    numData = length(model.t);
    %/~
    %   switch model.nonLinearity
    %    case 'linear'
    %     for k = intPoints
    %       for l = intPoints
    %         temp=0;
    %         for i = 1:length(model.t)
    %           arg1 = model.t(i)-model.mapt(k);
    %           arg2 = model.t(i)-model.mapt(l);
    %           if arg1 >= 0 & arg2 >= 0
    %             for j = 1:model.numGenes
    %               beta_ij = 1/model.yvar(i+(j-1)*length(model.t));
    %               temp = temp+beta_ij*exp(-model.D(j)*(arg1+arg2))*model.S(j)^2;
    %             end
    %           end
    %         end
    %         model.W(k, l) = temp;
    %       end
    %     end
    %     % Only needs to be done once for the linear system.
    %     model.updateW = false;
    %    otherwise
    %~/
    step2 = model.step*model.step;
    S2 = model.S.*model.S;
    for k = intPoints
        for l = intPoints
            temp=0;
            for i = 1:numData
                arg1 = model.t(i)-model.mapt(k);
                arg2 = model.t(i)-model.mapt(l);
                if arg1 >= 0 && arg2 >= 0
                    for j = 1:model.numGenes
                        if model.ngParam
                          gInd = j;
                        else
                          gInd = 1;
                        end
                        ind = i + (j-1)*numData;
                        beta_ij = 1/yvar(ind);
                        temp = temp+beta_ij*model.g_grad(k,gInd)* ...
                               model.g_grad(l,gInd)*exp(-model.D(j)* ...
                                                        (arg1+arg2)+ ...
                                                        log(S2(j))+log(step2));
                    end
                end
            end
            model.W(k, l) = temp;
            %/~
            if isinf(model.W(k, l)) | isnan(model.W(k, l))
                warning('model.W is inf/nan')
            end
            %~/
        end
    end
    for k=intPoints
        temp=0;
        for i=1:numData
            arg = model.t(i)-model.mapt(k);
            if arg >= 0
                for j=1:model.numGenes
                  if model.ngParam
                    gInd = j;
                  else
                    gInd = 1;
                  end
                  ind = i + (j-1)*numData;
                  beta_ij=1/yvar(ind);
                  factor=(model.ypred(model.times_index(i), j)- ...
                          model.y(ind))*beta_ij;
                  temp=temp + factor*model.g_grad2(k,gInd) * exp(- ...
                        model.D(j)*arg+log(model.S(j)) +log(model.step));
                end
            end
        end
        model.W(k, k) =  model.W(k, k) + temp;
        %/~
        if isinf(model.W(k, k)) | isnan(model.W(k, k))
            warning('model.W is inf/nan')
        end
        %~/
    end
    model.updateW = true;
    %  end
    model.W = 0.5*(model.W + model.W');
end
if strcmp(model.nonLinearity, 'linear')
    % Only needs to be done once for linear system.
    model.updateW = false;
end
