function [gK_uf, gK_ff] = gpsimCandidateCovGrads(model, M)

% GPSIMCANDIDATECOVGRADS Sparse objective function gradients wrt Covariance function.
% FORMAT
% DESC gives the gradients of the log likelihood with respect to the
% components of the posterior covariance.
% ARG model : the model for which the gradients are to be computed.
% ARG M : The training data for which the computation is to be made
% RETURN gK_uf : the gradient of the likelihood with respect to the
% elements of K_uf.
% RETURN gK_ff : the gradient of the likelihood with respect to K_ff
% 
% COPYRIGHT : Neil D. Lawrence, 2007
%
% SEEALSO : gpsimCreate, gpsimAddCandidate, gpsimCandidateLogLikeGradient

% SHEFFIELDML
%E = model.candidate.K_uf*M;
E = model.candidate.K_uf*model.candidate.invK*M;
EET = E*E';
AinvE = model.candidate.Ainv*E;
%AinvEET = model.candidate.Ainv*EET;
diagK_fuAinvEMT = sum(model.candidate.K_uf.*(model.candidate.Ainv*E*M'), 1)';
K_fuAinvEMT = model.candidate.K_uf'*model.candidate.Ainv*E*M';
AinvEETAinv = AinvE*AinvE';
K_ufdAinvplusAinvEETAinvK_fu = model.candidate.K_uf'*(model.candidate.Ainv ...
                                                  + AinvEETAinv)*model.candidate.K_uf;
invK_uuK_uf = model.invK*model.candidate.K_uf;
invK_uuK_ufDinv = invK_uuK_uf*model.candidate.invK;
MMT = M*M';
Q = -model.candidate.K + MMT ...
    + K_ufdAinvplusAinvEETAinvK_fu...
    -K_fuAinvEMT - K_fuAinvEMT';
%gK_uu = 0.5*((model.invK ...
%              -model.candidate.Ainv) - AinvEETAinv ...
%             + invK_uuK_ufDinv*Q*invK_uuK_ufDinv');
gK_uf = -invK_uuK_ufDinv*Q*model.candidate.invK ...      
        -model.candidate.Ainv*model.candidate.K_uf*model.candidate.invK ...
        -AinvEETAinv*model.candidate.K_uf*model.candidate.invK ...
        +model.candidate.Ainv*E*M'*model.candidate.invK;
gK_ff = 0.5*model.candidate.invK*Q*model.candidate.invK;
