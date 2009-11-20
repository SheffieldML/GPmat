function drosPreprocessExpData,

% DROSPREPROCESSEXPDATA Preprocess mmgMOS expression data further.
% FORMAT
% DESC Fits a normal distribution to absolute expression values
% from mmgMOS.
% ARG none : reads input from the current directory.
% RETURN none : writes output to the current directory.
%
% SEEALSO : drosLoadData
%
% COPYRIGHT : Antti Honkela, 2009

% DISIMRANK

FILEPREFIX = 'mmgmos_exprs';

expdata = importdata([FILEPREFIX, '_exprs.csv']);

normalisation = mean(expdata.data) - mean(mean(expdata.data));

prctiles = [5, 25, 50, 75, 95];
pcts = cell(5, 1);

for k=1:5,
  foo = importdata(sprintf('%s_prctile%d.csv', FILEPREFIX, prctiles(k)));
  pcts{k} = foo.data - repmat(normalisation, [size(foo.data, 1), 1]);
end
foo.rowheaders = foo.textdata(2:end, 1);

pctiles = cat(3, pcts{:});

fitmean = zeros(size(expdata.data));
fitvar = zeros(size(expdata.data));

for k=1:size(pctiles, 1),
  fprintf('%s\n', foo.rowheaders{k});
  prof = squeeze(pctiles(k, :, :));
  for l=1:36,
    t = distfit(exp(prof(l, :)), 'normal');
    fitmean(k, l) = t(1);
    fitvar(k, l) = t(2) .^ 2;
  end
end

fid = fopen(sprintf('%s_fitmean.txt', FILEPREFIX), 'w');
for k=1:size(pctiles, 1),
  fprintf(fid, '%s', foo.rowheaders{k});
  fprintf(fid, '\t%.14f', fitmean(k, :));
  fprintf(fid, '\n');
end
fclose(fid);
fid = fopen(sprintf('%s_fitvar.txt', FILEPREFIX), 'w');
for k=1:size(pctiles, 1),
  fprintf(fid, '%s', foo.rowheaders{k});
  fprintf(fid, '\t%.14f', fitvar(k, :));
  fprintf(fid, '\n');
end
fclose(fid);
