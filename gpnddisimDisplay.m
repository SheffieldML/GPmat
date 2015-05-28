function gpnddisimDisplay(model, spaceNum)

% GPNDDISIMDISPLAY Display a Gaussian process model.
% FORMAT
% DESC displays in human readable form the contents of the GPNDDISIM
% model.
% ARG model : the model structure to be displaced.
% ARG spaceNum : number of spaces to place before displaying model
% structure.
%
% SEEALSO : gpnddisimCreate, modelDisplay.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% COPYRIGHT : Antti Honkela, 2011


if nargin > 1
  spacing = repmat(32, 1, spaceNum);
else
  spaceNum = 0;
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Gaussian process non-decaying driving input single input motif model:\n')
fprintf(spacing);
fprintf('  Number of time points: %d\n', size(model.t, 1));
fprintf(spacing);
fprintf('  Number of genes: %d\n', model.numGenes);

fprintf(spacing);
fprintf('  SIM mean: %d\n', model.simMean);
fprintf(spacing);
fprintf('  DISIM start mean: %d\n', model.disimStartMean);

if any(model.mu~=0)
  fprintf(spacing);
  fprintf('  Output biases:\n');
  for i = 1:length(model.B)
    fprintf(spacing);
    fprintf('    Basal transcription gene %d: %2.4f\n', i, model.B(i));
  end
end

fprintf(spacing);
fprintf('  Kernel:\n')
kernDisplay(model.kern, 4+spaceNum);

