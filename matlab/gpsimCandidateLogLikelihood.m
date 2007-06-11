function ll = gpsimCandidateLogLikelihood(model)

% GPSIMCANDIDATELOGLIKELIHOOD Compute the log likelihood of a gene.
% FORMAT
% DESC computes the log likelihood of one gene given the existing genes
% in the GPSIM model. For use when ranking other prospective candidate genes.
% ARG model : the model containing the existing genes for which the log likelihood is computed.
% RETURN ll : the log likelihood of the new gene.
% 
% SEEALSO : gpsimCreate, gpsimAddCandidate, gpsimCandidateLogLikelihood,
% gpsimCandidateObjective
%
% COPYRIGHT : Neil D. Lawrence, 2007

% GPSIM


base = (model.logDetK - model.candidate.logDetA) + 
  
dim = size(model.y, 1);

invK_ffm = model.candiate.invK*model.m;
K_ufinvK_ffm = model.K_uf*invK_ffm;

ll = -0.5*(dim*log(2*pi) + model.candidate.logDetK ...
           + model.candidate.logDetA ...
           - model.logDetK ...
           + sum(sum(invK_ffm.*model.m))...
           - sum(sum((model.candidate.Ainv*K_ufinvK_ffm).*K_ufinvK_ffm)));

%/~ In case we need priors in.
%ll = ll + kernPriorLogProb(model.candidate.kern);
%ll = ll + priorLogProb(model.candidate.bprior, model.B);
%~/