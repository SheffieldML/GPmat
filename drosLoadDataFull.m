% DROSLOADDATA Loads the Drosophila data from the current directory
% FORMAT
% DESC Loads the Drosophila data from the current directory.
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if exist('drosophila_data.mat', 'file'),
  load('drosophila_data.mat');
elseif exist('data/drosophila_data.mat', 'file'),
  load('data/drosophila_data.mat');
else
  DATAPATH = './';

  expdata = importdata([DATAPATH, 'expro/mmgmos_exprs_exprs.csv']);
  normalisation = mean(expdata.data) - mean(mean(expdata.data));
  expdata.data = expdata.data - repmat(normalisation, [size(expdata.data, 1), 1]);

  expse = importdata([DATAPATH, 'expro/mmgmos_exprs_se.csv']);
  
  expfitm = importdata([DATAPATH, 'expro/mmgmos_exprs_fitmean.txt']);
  expfitv = importdata([DATAPATH, 'expro/mmgmos_exprs_fitvar.txt']);

  prctiles = [5, 25, 50, 75, 95];
  pcts = cell(5, 1);

  for k=1:5,
    foo = importdata(sprintf('%s/expro/mmgmos_exprs_prctile%d.csv', DATAPATH, prctiles(k)));
    pcts{k} = foo.data - repmat(normalisation, [size(foo.data, 1), 1]);
  end

  fid = fopen([DATAPATH, 'expro/dros_fbgn_annotations.txt']);
  FBgns = textscan(fid,'%s%s','Delimiter','\t');
  fclose(fid);

  fid = fopen([DATAPATH, 'expro/dros_symbol_annotations.txt']);
  Symbols = textscan(fid,'%s%s','Delimiter','\t');
  fclose(fid);

  drosexp.pctvalues = prctiles;
  drosexp.pctiles = cat(3, pcts{:});
  drosexp.probes = expdata.textdata(2:end, 1);
  drosexp.genes = FBgns{2};
  drosexp.symbols = Symbols{2};
  drosexp.geneids = drosMapGenesToIDs(drosexp.genes);
  
  drosexp.mean = expdata.data;
  drosexp.se = expse.data;
  drosexp.fitmean = expfitm.data;
  drosexp.fitvar = expfitv.data;

  chipdata = importdata(sprintf('%s/chipchip/chip_validation_distances.txt', DATAPATH));
  chipdistances.labels = chipdata.textdata(1, 2:end);
  chipdistances.data = chipdata.data;
  chipdistances.genes = chipdata.textdata(2:end, 1);
  chipdistances.geneids = drosMapGenesToIDs(chipdistances.genes);

  drosTF.names = {'tin', 'bin', 'twi', 'bap', 'mef2'};
  drosTF.probes = struct('tin', {'143426_at'}, 'bin', {'148296_at'}, ...
			 'twi', {'143396_at'}, 'bap', {'143523_at'}, ...
			 'mef2', {'153628_at'});
  drosTF.fbgns = struct('tin', {'FBgn0004110'}, 'bin', {'FBgn0045759'}, ...
			'twi', {'FBgn0003900'}, 'bap', {'FBgn0004862'}, ...
			'mef2', {'FBgn0011656'});

  insitudata = importdata([DATAPATH, 'insitu/insitu_data.txt']);
  drosinsitu.labels = insitudata.textdata(1, 2:end);
  drosinsitu.genes = insitudata.textdata(2:end, 1);
  drosinsitu.data = insitudata.data;
  drosinsitu.geneids = drosMapGenesToIDs(drosinsitu.genes);

  SIGNIFICANCE_LEVEL = .1;
  CHIPTHRESHOLD = 2000;

  tflabels = drosTF.names;
  
  drosmutant = [];

  for k=1:length(tflabels),
    tf = tflabels{k};

    I = drosFindGeneinds(chipdistances, drosexp.geneids, 1);
    J = strmatch(tf, chipdistances.labels);
    chip_validation.(tf) = NaN * ones(size(drosexp.genes));
    chip_validation.(tf)(I~=0) = chipdistances.data(I(I~=0), J) < CHIPTHRESHOLD;

    try,
      d = importdata(sprintf('%s/knockout/%s_knockout_qvalues.txt', DATAPATH, tflabels{k}));
      genes = d.rowheaders;
      drosmutant.(tf).genes = genes;
      drosmutant.(tf).geneids = drosMapGenesToIDs(genes);
      drosmutant.(tf).data = d.data < SIGNIFICANCE_LEVEL;
      drosmutant.(tf).qvalues = d.data;
    
      I = drosFindGeneinds(drosmutant.(tf), drosexp.geneids, 1);
      mutant_validation.(tf) = NaN * ones(size(drosexp.genes));
      mutant_validation.(tf)(I~=0) = drosmutant.(tf).data(I(I~=0));
    catch
    end
  end
  
  clear genes d expdata expse expfitm expfitv pcts prctiles foo chipdata insitudata FBgns Symbols
end
