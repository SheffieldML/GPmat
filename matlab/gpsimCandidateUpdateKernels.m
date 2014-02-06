function model = gpsimCandidateUpdateKernels(model)

% GPSIMCANDIDATEUPDATEKERNELS Updates the kernel representations in the GPSIM candidate structure.
% FORMAT
% DESC updates any representations of the kernel in the model's
% candidate sub-structure, such as the posterior covariance and its determinant.
% ARG model : the model structure for which candiate gene's kernels are being
% updated.
% RETURN model : the model structure with the candidate gene's kernels updated.
%
% SEEALSO gpsimCandidateExpandParam, gpsimCreate, gpsimUpdateKernels
%
% COPYRIGHT Neil D. Lawrence, 2007

% SHEFFIELDML


numBlocks = model.kern.numBlocks;
% for j = 1:model.candidate.numGenes
%   for i = 1:j-1
%     model.candiate.Kcomp{i, j} = 
%   end
%   model.candidate.Kcomp{j, j} = 
% end
model.candidate.K = zeros(size(model.candidate.y, 1));
numCanTimePoints = length(model.candidate.t);
numTimePoints = length(model.t);
iEndVal = 0;
for i = 1:model.candidate.numGenes
  iStartVal = iEndVal + 1;
  iEndVal = iEndVal + numCanTimePoints;
  % Compute diagonal components of the candidate points.
  model.candidate.K(iStartVal:iEndVal, iStartVal:iEndVal) = ...
      real(multiKernComputeBlock(model.candidate.kern, ...
                                 model.candidate.t, ...
                                 model.candidate.t, ...
                                 numBlocks + i, ...
                                 numBlocks + i));
  % Compute cross covariances with other candidate points.
  jEndVal = 0;
  for j = 1:i-1
    jStartVal = jEndVal + 1;
    jEndVal = jEndVal + numCanTimePoints;
    model.candidate.K(iStartVal:iEndVal, jStartVal:jEndVal) = ...
        real(multiKernComputeBlock(model.candidate.kern, ...
                                   model.candidate.t, ...
                                   model.candidate.t, ...
                                   numBlocks + i, ...
                                   numBlocks + j));
    model.candidate.K(jStartVal:jEndVal, iStartVal:iEndVal) = ...
        model.candidate.K(iStartVal:iEndVal, jStartVal:jEndVal)';
  end
  % Compute cross covariance with 'training set' genes.
  jEndVal = 0;
  for j = 1:model.numGenes
    jStartVal = jEndVal + 1;
    jEndVal = jEndVal + numTimePoints;
    model.candidate.K_uf(jStartVal:jEndVal, iStartVal:iEndVal) = ...
        real(multiKernComputeBlock(model.candidate.kern, ...
                                   model.t, ...
                                   model.candidate.t, ...
                                   numBlocks + i, ...
                                   j));
  end
end
model.candidate.K = model.candidate.K + diag(model.candidate.yvar);
[model.candidate.invK, U, jitter] = pdinv(model.candidate.K);
if jitter>1e-4
  fprintf('Warning: gpsimCandidateUpdateKernels added jitter of %2.4f to K\n', jitter)
end
model.candidate.logDetK = logdet(model.candidate.K, U);

% Update A
model.candidate.A = model.K ...
    + model.candidate.K_uf*model.candidate.invK*model.candidate.K_uf';

[model.candidate.Ainv, UA, jitter] = pdinv(model.candidate.A);
if jitter>1e-4
  fprintf('Warning: gpsimCandidateUpdate Kernels added jitter of %2.4f to A\n', jitter);
end
model.candidate.logDetA = logdet(model.candidate.A, UA);
