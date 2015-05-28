function genes = drosTargets(tf),

% DROSTARGETS returns a list of hand selected enchanced targets for given tf
%
% Usage:
%   genes = drosTargets(tf)
% where tf is one of 'bin', 'mef2', 'tin', 'twi'
%
% COPYRIGHT : Antti Honkela, 2007
  

% 'FBgn0029082', 'FBgn0036459', 
tingenes1 = {'FBgn0003117', 'FBgn0033459', 'FBgn0051368'};
tingenenames = {'pannier', 'CG12744', 'CG31368'};

tingenes2 = {'FBgn0003117', 'FBgn0015904', 'FBgn0030054', ...
	     'FBgn0027611', 'FBgn0029082'};

mefgenes1 = {'FBgn0000412', ...  % ???
	     'FBgn0000416', ...
	     'FBgn0000497', ...
	     'FBgn0001624', ...  % Twist bound nearby, probably for another gene
	     'FBgn0004401', ...
	     'FBgn0010109'};     % Lots of other activity nearby, probably for other genes

mefgenes2 = {'FBgn0031717', 'FBgn0035158', 'FBgn0036282', 'FBgn0004028', ... % !!!
	     'FBgn0032919', ... % This or CG1099
	     'FBgn0033688'};

twigenes1 = {'FBgn0003062', 'FBgn0003520', 'FBgn0014033', 'FBgn0025776', ...
	     'FBgn0028999', 'FBgn0029146'};
%twigenes1 = {'FBgn0003062', 'FBgn0003313', 'FBgn0004646', 'FBgn0003520', 'FBgn0014033', 'FBgn0020391'};

twigenes2 = {'FBgn0030749', ... % Mef2 and bin activity nearby
	     'FBgn0040108', ... % !!!
	     'FBgn0003062', ... % Mef2 nearby, probably for other genes
	     'FBgn0003308', 'FBgn0003607', 'FBgn0004635'}; % Good targets

bingenes1 = {'FBgn0004880', 'FBgn0010225', 'FBgn0015011', 'FBgn0026562', ...
	     'FBgn0033661', 'FBgn0035091', 'FBgn0038881', 'FBgn0039704'};

bingenes2 = {'FBgn0022073', 'FBgn0037240', 'FBgn0039268', 'FBgn0001197', ...
	     'FBgn0003890', 'FBgn0015011'};

bapgenes1 = {'FBgn0031349'};

switch tf,
 case 'twi',
  %genes = twigenes1([1 2 3 4 6]);
  genes = twigenes2;
 case 'mef2',
  %genes = mefgenes1([2 3 5]);
  %genes = mefgenes2([1 2 3 5 6]);
  genes = mefgenes2;
 case 'tin',
  %genes = tingenes1;
  genes = tingenes2;
 case 'bin',
  %genes = bingenes1(1:6);
  genes = bingenes2;
 %case 'bin1',
 % genes = bingenes1;
 otherwise,
  error('Unknown TF');
end
