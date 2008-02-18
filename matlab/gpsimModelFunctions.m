function gf = gpsimModelFunctions(model, options)
  
% GPSIMMODELFUNCTIONS Update the nonlinear transformation of f.
% FORMAT
% DESC updates the fields of model associated with the non-linear
% transformation of f, namely g, g_grad and g_grad2, etc.
% ARG model : the model to be updated.
% ARG model : the model with the updated g representation.
%
% COPYRIGHT : Pei Gao, Neil D. Lawrence, 2008
%
% SEEALSO : gpsimMapUpdateG

% GPSIM  
  
  
  if nargin < 2
    options = 'g';
  end

  if model.ngParam
    gf = zeros(length(model.f), model.numGenes);
    for k=1:model.numGenes
      gParamk = model.gParam(:,k);
      switch model.nonLinearity{k}
       case 'repression'
        gamma = gParamk(1);
        
        expf = expTransform(model.f, 'atox');
        exp2f = expf.*expf;
        exp3f = expf.*exp2f;        
        
        switch options
         case 'g'
          gf(:,k) = 1./(gamma+expf);
         case 'grad'
          gf(:,k) = -expf./((gamma+expf).*(gamma+expf));
         case 'grad2'
          grad(:,k) = -expf./((gamma+expf).*(gamma+expf));
          gf(:,k) = grad(:,k) + 2*exp2f./(gamma+expf).^3;
         case 'grad3'
          grad(:,k) = -expf./((gamma+expf).*(gamma+expf));
          grad2(:,k) = grad(:,k) + 2*exp2f./(gamma+expf).^3;
          gf(:,k) = grad2(:,k) + 4*exp2f./(gamma+expf).^3 - ...
                    6*exp3f./(gamma+expf).^4;
         case 'paramGrad'
          gf(:,k) = -1./(gamma+expf).^2;
         case 'paramGgrad'
          gf(:,k) = 2*expf./((gamma+expf).^3);
         case 'paramGgrad2'
          gf(:,k) = 2*expf./((gamma+expf).^3) - 6*exp2f./((gamma+ ...
                                                            expf).^4);
        end    


%        case 'repression'
%         gamma = gParamk(1);
%         
%         f = model.f;
%       
%         switch options
%          case 'g'
%           gf(:,k) = 1./(gamma+f.*f);
%          case 'grad'
%           gf(:,k) = -2*f./((gamma+f.*f).*(gamma+f.*f));
%          case 'grad2'
%           gf(:,k) = -2./((gamma+f.*f).*(gamma+f.*f))+8*f.*f./(gamma+f.*f).^3;
%          case 'grad3'
%           gf(:,k) = 24*f./(gamma+f.*f).^3 - 48*f.^3./(gamma+f.*f).^4;
%          case 'paramGrad'
%           gf(:,k) =  -1./(gamma+f.*f).^2;
%          case 'paramGgrad'
%           gf(:,k) = 4*f./((gamma+f.*f).^3);
%          case 'paramGgrad2'
%           gf(:,k) = 4*f./((gamma+f.*f).^3) - 24*f.*f./((gamma+f.*f).^4);
%         end                   
        
  
       case 'activation'
        gamma = gParamk(1);
         
        expmf = expTransform(-model.f, 'atox');
        expm2f = expmf.*expmf;
        expm3f = expmf.*expm2f;
        switch options
         case 'g'
          gf(:,k) = 1./(1+gamma*expmf);
         case 'grad'
          gf(:,k) = gamma*expmf./((1+gamma*expmf).*(1+gamma*expmf));
         case 'grad2'
          grad(:,k) = -gamma*expmf./((1+gamma*expmf).*(1+gamma*expmf));
          gf(:,k) = grad(:,k) + 2*gamma^2*expm2f./(1+gamma*expmf).^3;
         case 'grad3'
          grad(:,k) = -gamma*expmf./((1+gamma*expmf).*(1+gamma*expmf));
          grad2(:,k) = grad(:,k) + 2*gamma^2*expm2f./(1+gamma*expmf).^3;
          gf(:,k) = -grad2(:,k) - 4*gamma^2*expm2f./(1+gamma*expmf).^3 ...
                    + 6*gamma^3*expm3f./(1+gamma*expmf).^4;
         case 'paramGrad'
          gf(:,k) = -expmf./((1+gamma*expmf).*(1+gamma*expmf));
         case 'paramGgrad'
          grad = expmf./((1+gamma*expmf).^2);
          gf(:,k) = grad - 2*gamma*expm2f./((1+gamma*expmf).^3);
         case 'paramGgrad2'
          gf(:,k) = expmf./((1+gamma*expmf).^2);
          gf(:,k) = gf(:,k) - 2*gamma*expm2f./((1+gamma*expmf).^3);
          gf(:,k) = -gf(:,k) + 4*gamma*expm2f./((1+gamma*expmf).^3) ...
                    - 6*gamma^2*expm3f./((1+gamma*expmf).^4);
        end    
      end
    end
  end
      
    
  
