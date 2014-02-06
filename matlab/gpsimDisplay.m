function gpsimDisplay(model, spaceNum)

% GPSIMDISPLAY Display a Gaussian process model.
% FORMAT
% DESC displays in human readable form the contents of the GPSIM
% model.
% ARG model : the model structure to be displaced.
% ARG spaceNum : number of spaces to place before displaying model
% structure.
%
% SEEALSO : gpsimCreate, modelDisplay.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% SHEFFIELDML

if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Gaussian process single input motif model:\n')
fprintf(spacing);
fprintf('  Number of time points: %d\n', size(model.t, 1));
fprintf(spacing);
fprintf('  Number of genes: %d\n', model.numGenes);

if any(model.mu~=0)
  fprintf(spacing);
  fprintf('  Output biases:\n');
  for i = 1:length(model.mu)
    fprintf(spacing);
    fprintf('    Basal transcription gene %d: %2.4f\n', i, model.B(i));
  end
end

fprintf(spacing);
fprintf('  Kernel:\n')
kernDisplay(model.kern, 4+spaceNum);

