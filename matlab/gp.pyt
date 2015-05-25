def [gmu, gsigmavar] = gpPosteriorGradMeanVar(model, X):

	"""Gadient of the mean and variances of the posterior at points given by X.
	
	Description:
	
	[gmu, gsigmavar] = gpPosteriorGradMeanVar(model, X) computes the
	 gradient of the mean and variances of the posterior distribution
	 of a Gaussian process with respect to the input locations.
	 Returns:
	  gmu - the gradient of the posterior mean with respect to the input
	   locations.
	  gsigmavar - the gradient of the posterior variances with respect
	   to the input locations.
	 Arguments:
	  model - the model for which gradients are to be computed.
	  X - the input locations where gradients are to be computed.
		

	See also
	GPCREATE, GPPOSTERIORMEANVAR


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	if ~isfield(model, 'alpha')
	  model = gpComputeAlpha(model);
	end
	
	if size(X, 1) > 1
	  error('This function only handles one data-point at a time')
	end
	
	switch model.approx
	 case 'ftc'
	  gX = kernGradX(model.kern, X, model.X);
	  kX = kernCompute(model.kern, X, model.X)';
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  gX = kernGradX(model.kern, X, model.X_u);
	  kX = kernCompute(model.kern, X, model.X_u)';
	 otherwise
	  error('Unrecognised approximation type');
	end
	
	diaggK = kernDiagGradX(model.kern, X);
	
	gmu = zeros(size(X, 2), model.d);
	gsigmavar = zeros(size(X, 2), model.d);
	
	switch model.approx
	 case 'ftc'
	  Kinvgk = model.invK_uu*gX;
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  Kinvgk = (model.invK_uu - (1/model.beta)*model.Ainv)*gX;
	 otherwise
	  error('Unrecognised approximation type');
	end
	
	gsigmavar = repmat(diaggK' - 2*Kinvgk'*kX, 1, model.d);
	gmu = gX'*model.alpha;
	
	gmu = gmu.*repmat(model.scale, model.q, 1);
	gsigmavar = gsigmavar.*repmat(model.scale.*model.scale, model.q, 1);
def model = gpCovGradsTest(model):

	"""Test the gradients of the likelihood wrt the covariance.
	
	Description:
	
	model = gpCovGradsTest(model) tests the gradients of the
	 covariance to ensure they are correct.
	 Returns:
	  model - the model that was tested.
	 Arguments:
	  model - the model to be tested.
		

	See also
	GPCREATE, GPCOVGRADS


	Copyright (c) 2006, 2009 Neil D. Lawrence
	
	"""
		
	% WARNING --- this isn't testing g_Lambda in gpCovGrads
	  
	changeVal = 1e-6;
	switch model.approx
	 case 'ftc'
	  
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  for i =1 :size(model.K_uu, 1)
	    for j=1:i
	      origK = model.K_uu(i, j);
	      model.K_uu(i, j) = origK + changeVal;
	      model.K_uu(j, i) = model.K_uu(i, j);
	      [model.invK_uu, U] = pdinv(model.K_uu);
	      model.logDetK_uu = logdet(model.K_uu, U);
	      model = gpUpdateAD(model);
	      objPlus = gpLogLikelihood(model);
	      model.K_uu(i, j) = origK - changeVal;
	      model.K_uu(j, i) = model.K_uu(i, j);
	      [model.invK_uu, U] = pdinv(model.K_uu);
	      model.logDetK_uu = logdet(model.K_uu, U);
	      model = gpUpdateAD(model);
	      objMinus = gpLogLikelihood(model);
	      diffsK_uu(i, j) = (objPlus - objMinus)/(2*changeVal);
	      diffsK_uu(j, i) = diffsK_uu(i, j);
	      model.K_uu(i, j) = origK;
	      model.K_uu(j, i) = origK;
	      [model.invK_uu, U] = pdinv(model.K_uu);
	      model.logDetK_uu = logdet(model.K_uu, U);
	      model = gpUpdateAD(model);
	    end
	  end
	  for i=1:size(model.K_uf, 1)
	    for j=1:size(model.K_uf, 2)
	      origK = model.K_uf(i, j);
	      model.K_uf(i, j) = origK + changeVal;
	      model = gpUpdateAD(model);
	      objPlus = gpLogLikelihood(model);
	      model.K_uf(i, j) = origK - changeVal;
	      model = gpUpdateAD(model);
	      objMinus = gpLogLikelihood(model);
	      diffsK_uf(i, j) = (objPlus - objMinus)/(2*changeVal);
	      model.K_uf(i, j) = origK;
	      model = gpUpdateAD(model);
	    end
	  end
	  
	  [gK_uu, gK_uf, g_Lambda] = gpCovGrads(model, model.m);
	  
	  gK_uuMaxDiff = max(max(abs(2*(gK_uu-diag(diag(gK_uu))) ...
	                             + diag(diag(gK_uu)) ...
	                             - diffsK_uu)));
	  gK_ufMaxDiff = max(max(abs(gK_uf - diffsK_uf)));
	  
	  fprintf('K_uu grad max diff %2.4f\n', gK_uuMaxDiff);
	  if gK_uuMaxDiff > 1e-4
	    disp(2*(gK_uu-diag(diag(gK_uu))) ...
	         + diag(diag(gK_uu)) ...
	         - diffsK_uu);
	  end
	  fprintf('K_uf grad max diff %2.4f\n', gK_ufMaxDiff);
	  if gK_ufMaxDiff > 1e-4
	    disp(gK_uf - diffsK_uf)
	  end
	end
	

def y = gpSubspaceOut(model,x):

	""
	
	Description:
		


	Copyright (c) 2008 Carl Henrik Ek
	
	"""
		  
	y = NaN.*ones(size(x,1),length(model.dim));
	y(:,find(model.dim)) = gpOut(model,x);
	
	return;
def model = gpUpdateKernels(model, X, X_u):

	"""Update the kernels that are needed.
	
	Description:
	
	model = gpUpdateKernels(model) updates any representations of the
	 kernel in the model structure, such as invK, logDetK or K.
	 Returns:
	  model - the model structure with the kernels updated.
	 Arguments:
	  model - the model structure for which kernels are being updated.
		DESC updates any representations of the kernel in the model
		structure, such as invK, logDetK or K.
		ARG model : the model structure for which kernels are being
		updated.
		ARG X : the input locations for update of kernels.
		ARG X_u : the inducing input locations.
		RETURN model : the model structure with the kernels updated.
		
		

	See also
	GPEXPANDPARAM, GPCREATE


	Copyright (c) 2005, 2006, 2007, 2009 Neil D. Lawrence
	
	"""
		
	jitter = 1e-6;
	
	switch model.approx
	 case 'ftc'
	  % Long term should allow different kernels in each dimension here.
	  model.K_uu = kernCompute(model.kern, X);
	    
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    % Add inverse beta to diagonal if it exists.
	    if isfield(model, 'beta') & ~isempty(model.beta)
	      model.K_uu(1:size(model.K_uu, 1)+1:end) = ...
	          model.K_uu(1:size(model.K_uu, 1)+1:end) + 1./model.beta';
	    end
	    [model.invK_uu, U] = pdinv(model.K_uu);
	    model.logDetK_uu = logdet(model.K_uu, U);
	  else   
	    for i = 1:model.d
	      if isfield(model, 'beta') & ~isempty(model.beta)
	        if size(model.beta, 2) == model.d
	          betaAdd = model.beta(:, i)';
	        else
	          betaAdd = model.beta;
	        end
	        if isfield(model, 'beta') & ~isempty(model.beta)
	          model.K_uu(1:size(model.K_uu, 1)+1:end) = ...
	              model.K_uu(1:size(model.K_uu, 1)+1:end) + 1./betaAdd;
	        end
	      end
	      ind = gpDataIndices(model, i);
	      [model.invK_uu{i}, U] = pdinv(model.K_uu(ind, ind));
	      model.logDetK_uu(i) = logdet(model.K_uu(ind, ind), U);
	    end
	  end
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  model.K_uu = kernCompute(model.kern, X_u);
	  
	  if ~isfield(model.kern, 'whiteVariance') | model.kern.whiteVariance == 0
	    % There is no white noise term so add some jitter.
	    model.K_uu = model.K_uu ...
	        + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
	  end
	  model.K_uf = kernCompute(model.kern, X_u, X);
	  [model.invK_uu, model.sqrtK_uu] = pdinv(model.K_uu);
	  model.logDetK_uu = logdet(model.K_uu, model.sqrtK_uu);
	
	end
	
	switch model.approx
	 case {'dtcvar', 'fitc'}
	  model.diagK = kernDiagCompute(model.kern, X);
	 case {'pitc'}
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      model.K{i} = kernCompute(model.kern, X(ind, :));
	    end
	  else
	    for j = 1:model.d
	      for i = 1:length(model.blockEnd)
	        ind = gpDataIndices(model, j, i);
	        model.K{i, j} = kernCompute(model.kern, X(ind, :));
	      end
	    end
	  end
	
	end
	
	model = gpUpdateAD(model, X);  

def ll = gpPointLogLikelihood(model, x, y):

	"""Log-likelihood of a test point for a GP.
	
	Description:
	
	gpPointLogLikelihood(model, x, y) returns the log likelihood of a
	 latent point and an observed data point for the posterior
	 prediction of a GP model.
	 Arguments:
	  model - the model for which the point prediction will be made.
	  x - the input point for which the posterior distribution will be
	   evaluated.
	  y - the target point for which the posterior distribution will be
	   evaluated.
		

	See also
	GPCREATE


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	logTwoPi = log(2*pi);
	[mu, varSigma] = gpPosteriorMeanVar(model, x);
	ll = zeros(size(x, 1), 1);
	ydiff = y-mu;
	ll = log(varSigma) + (ydiff.*ydiff)./varSigma +logTwoPi;
	ll(find(isnan(ll)))=0;
	ll = -0.5*sum(ll, 2);


	"""Simple demonstration of sampling from a covariance function.
	
	Description:
	
	"""
		
	randn('seed', 1e5)
	rand('seed', 1e5)
	
	bw = false;
	x = linspace(-1, 1, 25)';
	kern = kernCreate(x, 'rbf');
	kern.inverseWidth = 10;
	
	figure(1)
	clf
	K = kernCompute(kern, x);
	imagesc(K);
	if bw
	  colormap gray
	else
	  colormap jet
	end  
	colorbar
	t = [];
	t = [t xlabel('n')];
	t = [t ylabel('m')];
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 24)
	if exist('printDiagram', 'var') && printDiagram
	  printPlot('gpCovariance', '../tex/diagrams/', '../html');
	end
	% need to take the real part of the sample as the kernel is numerically less than full rank 
	f = real(gsamp(zeros(1, size(x, 1)), K, 1))';
	
	figure(2) 
	clf
	a = plot(f, 'k.');
	t = [t text(12.5, -.5, 'n')];
	t = [t text(-5, .5, 'f_n')];
	set(t, 'fontname', 'times')
	set(t, 'fontsize', 24)
	set(t, 'fontangle', 'italic')
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 18)
	set(a,'markersize', 20)
	set(a, 'linewidth', 2)
	zeroAxes(gca)
	if exist('printDiagram', 'var') && printDiagram
	  printPlot('gpSample', '../tex/diagrams/', '../html');
	end
	% if bw
	%   set(a, 'color', [0 0 0])
	%   print('-deps', '../tex/diagrams/gpSamplebw.eps');
	% else
	%   print('-depsc', '../tex/diagrams/gpSample.eps');
	% end
	
	save demGpSample K x f
def [analMu, analCov, diffMu, diffCov] = gpPosteriorMeanCovarTest(model, X):

	"""Test the gradients of the mean and covariance.
	
	Description:
	
	[analMu, analCov, diffMu, diffCov, delta] =
	 gpPosteriorMeanCovarTest(model, X) tests the gpPosteriorMeanCovar
	 and gpPosteriorGradMeanCovar functions.
	 Returns:
	  analMu - the analytical gradients of the mean with respect to X.
	  analCov - the analytical gradients of the covariance with respect
	   to X.
	  diffMu - the numerical gradients of the mean with respect to X.
	  diffCov - the numerical gradients of the covariance with respect
	   to X.
	  delta - erros between the numerical and analytical gradients.
	 Arguments:
	  model - the model to test the gradients for.
	  X - the input locations to test the gradients for.
		

	See also
	GPPOSTERIORMEANCOVAR, GPPOSTERIORGRADMEANCOVAR


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	[analMu, analCov] = gpPosteriorGradMeanCovar(model, X);
	origX = X;
	change = 1e-6;
	for i = 1:size(X, 1)
	  for j = 1:size(X, 2)
	    X(i, j) = origX(i, j) - change;
	    [muMinus, covMinus] = gpPosteriorMeanCovar(model, X);
	    X(i, j) = origX(i, j) + change;
	    [muPlus, covPlus] = gpPosteriorMeanCovar(model, X);
	    X(i, j) = origX(i, j);
	    diffMu{j}(i, :) = (muPlus(i, :) - muMinus(i, :))/(2*change);
	    for k = 1:model.d
	      % Not sure why the transpose is need here ... hope it isn't
	      % two wrongs making a right ...
	      diffCov{j, k}(:, i) = (covPlus{k}(i, :) - covMinus{k}(i, :))'/(2*change);
	    end
	    
	  end
	end

	"""Do a simple 1-D regression after Snelson & Ghahramani's example.
	
	Description:
	
	"""
		
	% Fix seeds
	randn('seed', 2e5);
	rand('seed', 2e5);
	seedVal = 2e5;
	dataSetName = 'spgp1d';
	experimentNo = 2;
	
	% load data
	[X, y] = mapLoadData(dataSetName, seedVal);
	
	% Set up model
	options = gpOptions('fitc');
	options.numActive = 9;
	options.optimiser = 'conjgrad';
	
	% use the deterministic training conditional.
	q = size(X, 2);
	d = size(y, 2);
	
	model = gpCreate(q, d, X, y, options);
	model.X_u = randn(9, 1)*0.25 - 0.75;
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);
	
	% Optimise the model.
	iters = 1000;
	display = 1;
	model.beta = 4/var(y);
	model.kern.variance = var(y);
	model.kern.inverseWidth = 1./((-min(X)+max(X))'/2).^2
	
	model = gpOptimise(model, display, iters);
	
	% Save the results.
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	save(['dem' capName num2str(experimentNo) '.mat'], 'model');
	
	
	demSpgp1dPlot


	"""Load in the relevant toolboxes for GP.
	
	Description:
	
	"""
		importLatest('netlab');
	importLatest('mocap');
	importLatest('ndlutil');
	importLatest('prior');
	importLatest('mltools');
	importLatest('optimi');
	importLatest('datasets');
	importLatest('kern');
	importLatest('noise'); % Only needed for C++ load ins.
def demGpCovFuncSample

	"""Sample from some different covariance functions.
	
	Description:
	
	demGpCovFuncSample samples from some different covariance
	 functions with different kernel parameters.
		

	See also
	KERNCREATE


	Copyright (c) 2006, 2008 Neil D. Lawrence
	
	"""
		
	  global printDiagram
	  randn('seed', 1e5)
	  rand('seed', 1e5)
	  
	  x = linspace(-2, 2, 200)';
	  kern = kernCreate(x, 'rbf');
	  
	  numSamps = 10;
	  
	  figure(1);
	  kern.inverseWidth = 1/(0.3^2);
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	  
	  figure(2)
	  kern.inverseWidth = 1;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	  
	  figure(3)
	  kern.inverseWidth = 1/(0.3^2);
	  kern.variance = 4;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	  
	  figure(4)
	  kern = kernCreate(x, 'bias');
	  kern.variance = 4;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	
	  
	  figure(5)
	  kern = kernCreate(x, 'lin');
	  kern.variance = 16;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	
	  figure(6)
	  kern = kernCreate(x, 'poly');
	  kern.variance = 1;
	  kern.biasVariance = 1;
	  kern.weightVariance = 1;
	  kern.degree = 5;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	
	  figure(7)
	  kern = kernCreate(x, 'mlp');
	  kern.variance = sqrt(pi)/2;
	  kern.biasVariance = 100;
	  kern.weightVariance = 100;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	  
	  
	  figure(8)
	  kern = kernCreate(x, 'mlp');
	  kern.variance = sqrt(pi/2);
	  kern.biasVariance = 0;
	  kern.weightVariance = 100;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	  
	  
	  figure(9)
	  kern = kernCreate(x, {'rbf', 'bias', 'white'});
	  kern.comp{1}.variance = 1;
	  kern.comp{1}.inverseWidth = 10;
	  kern.comp{2}.variance = 1;
	  kern.comp{3}.variance = 0.01;
	  K = kernCompute(kern, x);
	  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
	  a = plot(x, f);
	  prepPlot(a);
	
	end
	
	
def prepPlot(a):
	
	  global printDiagram
	  set(gca, 'xlim', [-2 2])
	  set(gca, 'ylim', [-4 4])
	  set(gca, 'fontname', 'times')
	  set(gca, 'fontsize', 18)
	  set(a,'markersize', 10)
	  set(a, 'linewidth', 2)
	  zeroAxes(gca, 0.025, 18, 'times');
	  fig = gcf;
	  if printDiagram
	    printPlot(['demGpCovFuncSample' num2str(fig)], '../tex/diagrams', '../html');
	  end
	end
def f = gpObjective(params, model):

	"""Wrapper function for GP objective.
	
	Description:
	
	f = gpObjective(params, model) returns the negative log likelihood
	 of a Gaussian process model given the model structure and a vector
	 of parameters. This allows the use of NETLAB minimisation
	 functions to find the model parameters.
	 Returns:
	  f - the negative log likelihood of the GP model.
	 Arguments:
	  params - the parameters of the model for which the objective will
	   be evaluated.
	  model - the model structure for which the objective will be
	   evaluated.
		

	See also
	SCG, CONJGRAD, GPCREATE, GPGRADIENT, GPLOGLIKELIHOOD, GPOPTIMISE


	Copyright (c) 2005, 2006 Neil D. Lawrence
	
	"""
		
	model = gpExpandParam(model, params);
	f = - gpLogLikelihood(model);


	"""Do a simple 1-D regression after Snelson & Ghahramani's example.
	
	Description:
	
	"""
		
	% Fix seeds
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	dataSetName = 'spgp1d';
	experimentNo = 5;
	
	% load data
	[X, y] = mapLoadData(dataSetName);
	
	% Set up model
	options = gpOptions('dtcvar');
	options.kern = {'rbf', 'bias', 'white'}
	options.numActive = 9;
	options.optimiser = 'scg';
	
	% use the deterministic training conditional.
	q = size(X, 2);
	d = size(y, 2);
	
	model = gpCreate(q, d, X, y, options);
	model.kern.comp{3}.setVariance(1e-4);
	model.X_u = randn(options.numActive, 1)*0.1;
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);
	
	% Optimise the model.
	iters = 1000;
	display = 1;
	
	model = gpOptimise(model, display, iters);
	
	% Save the results.
	capName = dataSetName;
	capName(1) = upper(capName(1));
	save(['dem' capName num2str(experimentNo) '.mat'], 'model');
	
	
	demSpgp1dPlot


	"""Model silhouette data with independent RBF GPs.
	
	Description:
	
	"""
		% FORMAT
	% DESC runs a simple regression on the Agawal and Triggs data.
	%
	% SEEALSO : demSilhouetteGp1, demSilhouetteAverage
	% 
	% COPYRIGHT : Neil D. Lawrence, 2008
	
	
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	dataSetName = 'silhouette';
	experimentNo = 1;
	
	% load data
	[X, y, XTest, yTest] = mapLoadData(dataSetName);
	
	% Set up the model
	options = gpOptions('ftc');
	
	% Scale outputs to variance 1.
	options.scale2var1 = true;
	
	% Use the full Gaussian process model.
	q = size(X, 2);
	d = size(y, 2);
	model = gpCreate(q, d, X, y, options);
	
	display = 1;
	iters = 1000;
	
	model = gpOptimise(model, display, iters);
	modelDisplay(model)
	
	% Save results
	capName = dataSetName;
	capName(1) = upper(capName(1));
	fileBaseName = ['dem' capName 'Gp' num2str(experimentNo)];
	save([fileBaseName '.mat'], 'model');
	demSilhouettePlot

def m = gpComputeM(model):

	"""Compute the matrix m given the model.
	
	Description:
	
	gpComputeM(model, m) computes the matrix m (the scaled, bias and
	 mean function removed matrix of the targets), given the model.
	 Arguments:
	  model - the model for which the values are to be computed.
	  m - the scaled, bias and mean function removed values.
		

	See also
	GPCREATE, GPCOMPUTEALPHA, GPUPDATEAD


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	% Remove mean function value from m (if mean function present).
	if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
	  m = model.y - modelOut(model.meanFunction, model.X);
	else
	  m = model.y;
	end
	
	% Remove bias and apply scale.
	for i = 1:model.d
	  m(:, i) = m(:, i) - model.bias(i);
	  if model.scale(i)
	    m(:, i) = m(:, i)/model.scale(i);
	  end
	end

def model = gpComputeAlpha(model, m):

	"""Update the vector `alpha' for computing posterior mean quickly.
	
	Description:
	
	model = gpComputeAlpha(model, m) updates the vectors that are
	 known as `alpha' in the support vector machine, in other words
	 invK*y, where y is the target values.
	 Returns:
	  model - the model with the updated alphas.
	 Arguments:
	  model - the model for which the alphas are going to be updated.
	  m - the values of m for which the updates will be made.
		

	See also
	GPCREATE, GPUPDATEAD, GPUPDATEKERNELS


	Copyright (c) 2006, 2009 Neil D. Lawrence
	
	"""
		
	if nargin < 2
	  m = model.m;
	end
	
	switch model.approx
	 case 'ftc'
	  model.alpha = zeros(model.N, model.d);
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    model.alpha = model.invK_uu*m;
	  else
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      model.alpha(ind, i) = model.invK_uu{i}* ...
	          m(ind, i);
	    end
	  end
	 
	 case {'dtc', 'dtcvar'}
	  model.alpha = zeros(model.k, model.d);
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    model.alpha = model.Ainv*model.K_uf*m;
	  else
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      model.alpha(:, i) = model.Ainv{i} ...
	          *model.K_uf(:, ind) ...
	          *m(ind, i);
	    end
	  end
	 case 'fitc'
	  model.alpha = zeros(model.k, model.d);
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    model.alpha = model.Ainv*model.K_uf*model.Dinv*m;
	  else
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      model.alpha(:, i) = model.Ainv{i} ...
	          *model.K_uf(:, ind) ...
	          *model.Dinv{i}*m(ind, i);
	    end
	  end
	 case 'pitc'
	  model.alpha = zeros(model.k, model.d);
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      model.alpha = model.alpha+model.Ainv*model.K_uf(:, ind)* ...
	          model.Dinv{i}*m(ind, :);
	    end
	  else
	    for i = 1:length(model.blockEnd)
	      for j = 1:model.d
	        ind = gpDataIndices(model, j, i);
	        model.alpha(:, j) = model.alpha(:, j)+model.Ainv{j}*model.K_uf(:, ind)* ...
	            model.Dinv{i, j}*m(ind, j);
	      end
	    end  
	  end
	end
def model = gpReadFromFile(fileName, varargin):

	"""Load a file produced by the C++ implementation.
	
	Description:
	
	model = gpReadFromFile(fileName) loads a GP model from a file
	 produced by the C++ GP implementation.
	 Returns:
	  model - a MATLAB GP model structure containing the model from the
	   file.
	 Arguments:
	  fileName - the file name written by the C++ software.
		

	See also
	GPREADFROMFID, GPCREATE


	Copyright (c) 2005 Neil D. Lawrence
	
	"""
		
	FID = fopen(fileName);
	if FID==-1
	  error(['Cannot find file ' fileName])
	end
	model = gpReadFromFID(FID, varargin{:});
	fclose(FID);

	"""Demonstrate Gaussian processes for interpolation.
	
	Description:
	
	demInterpolation runs a simple one-D Gaussian process displaying
	 errorbars.
		

	See also
	GPCREATE, DEMREGRESSION


	Copyright (c) 2006, 2008 Neil D. Lawrence
	
	"""
		
	randn('seed', 1e6)
	rand('seed', 1e6)
	
	% Create data set
	x = linspace(-1, 1, 9)';
	trueKern = kernCreate(x, 'rbf');
	K = kernCompute(trueKern, x);
	% Sample some true function values.
	yTrue = gsamp(zeros(size(x))', K, 1)';
	
	markerSize = 20;
	markerWidth = 6;
	markerType = 'k.';
	lineWidth = 2;
	% Create a test set
	indTrain{1} = [1 9]';
	indTrain{2} = [1 5 9]';
	indTrain{3} = [1 3 5 7 9]';
	indTrain{4} = [1 2 3 4 5 6 7 8 9]';
	figNo = 1;
	fillColor = [0.7 0.7 0.7];
	for i = 0:length(indTrain)
	  if i > 0
	    yTrain = yTrue(indTrain{i});
	    xTrain = x(indTrain{i});
	    kern = kernCreate(x, 'rbf');
	    % Change inverse variance (1/(lengthScale^2)))
	    kern.inverseWidth = 5;
	    
	    xTest = linspace(-2, 2, 200)';
	    
	    Kx = kernCompute(kern, xTest, xTrain);
	    Ktrain = kernCompute(kern, xTrain, xTrain);
	    
	    yPred = Kx*pdinv(Ktrain)*yTrain;
	    yVar = kernDiagCompute(kern, xTest) - sum(Kx*pdinv(Ktrain).*Kx, 2);
	    ySd = sqrt(yVar);
	    figure(figNo)
	    clf
	    fill([xTest; xTest(end:-1:1)], ...
	         [yPred; yPred(end:-1:1)] ...
	         + 2*[ySd; -ySd], ...
	         fillColor,'EdgeColor',fillColor)
	    hold on;
	    h = plot(xTest, yPred, 'k-');
	    set(h, 'linewidth', lineWidth)
	    p = plot(xTrain, yTrain, markerType);
	    set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	    set(gca, 'xtick', [-2 -1 0 1 2]);
	    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
	    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
	    zeroAxes(gca);
	    if exist('printDiagram') && printDiagram
	      printPlot(['demInterpolation' num2str(figNo)], '../tex/diagrams', '../html');
	    end
	    figNo = figNo + 1;
	  else
	    p = [];
	  end
	  if i < length(indTrain)
	    figure(figNo)
	    if i>0
	      fill([xTest; xTest(end:-1:1)], ...
	           [yPred; yPred(end:-1:1)] ...
	           + 2*[ySd; -ySd], ...
	           fillColor,'EdgeColor',fillColor)
	      hold on
	      h = plot(xTest, yPred, 'k-');
	      set(h, 'linewidth', lineWidth)
	    end
	    p = [p plot(x(indTrain{i+1}), yTrue(indTrain{i+1}), markerType)];
	    set(p, 'markersize', markerSize, 'linewidth', markerWidth);
	    set(gca, 'xtick', [-2 -1 0 1 2]);
	    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
	    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
	    zeroAxes(gca);
	    if exist('printDiagram') && printDiagram
	      printPlot(['demInterpolation' num2str(figNo)], '../tex/diagrams', '../html');
	    end
	    figNo = figNo + 1;
	  end
	end
def g = gpScaleBiasGradient(model):

	"""Compute the log likelihood gradient wrt the scales.
	
	Description:
	
	gpScaleBiasGradient(model, g) computes the gradient of the log
	 likelihood with respect to the scales. In the future the gradients
	 with respect to the biases may also be included.
	 Arguments:
	  model - the model for which the gradients are to be computed.
	  g - the gradients of the likelihood with respect to the scales.
		

	See also
	GPCREATE, GPLOGLIKEGRADIENTS, GPLOGLIKELIHOOD


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	g = [];
	if model.learnScales
	  g = 1./model.scale.*(model.innerProducts-1);
	  fhandle = str2func([model.scaleTransform 'Transform']);
	  g = g.*fhandle(model.scale, 'gradfact');
	end


	""
	
	Description:
	
	"""
		
	% Show prediction for test data.
	yPred = modelOut(model, XTest);
	xyzankurAnimCompare(yPred, yTest);
	
	yDiff = (yPred - yTest);
	rmsError = sqrt(sum(sum(yDiff.*yDiff))/prod(size(yDiff)));
	
	counter = 0;
	if printDiagram
	  ind = 1:27:size(yPred, 1)
	  for i = ind
	    counter = counter + 1;
	    figure
	    handle = xyzankurVisualise(yPred(i,:), 1);
	    printPlot([fileBaseName '_' num2str(counter)], '../tex/diagrams', '../html') 
	  end
	end

def params = gpSubspaceExtractParam(model):

	""
	
	Description:
		


	Copyright (c) 2008 Carl Henrik Ek
	
	"""
		  
	params = gpExtractParam(model);
	
	return;

	"""Plot results from 1-D sparse GP.
	
	Description:
	
	"""
		
	fillColor = [0.7 0.7 0.7];
	xTest = linspace(-1.5, 1.5, 200)';
	[mu, varSigma] = gpPosteriorMeanVar(model, xTest);
	
	figure
	fill([xTest; xTest(end:-1:1)], ...
	     [mu; mu(end:-1:1)] ...
	     + 2*[sqrt(varSigma); -sqrt(varSigma(end:-1:1))], ...
	     fillColor,'EdgeColor',fillColor)
	hold on;
	plot(X, y, 'k.');
	a = plot(xTest, mu, 'k-');
	if isfield(model, 'X_u') && ~isempty(model.X_u)
	  b = plot(model.X_u, -ones(size(model.X_u)), 'bx');
	  set(b, 'linewidth', 2)
	  set(b, 'markersize', 10);
	end
	set(gca, 'ylim', [-1 2])
	set(gca, 'xlim', [-1.5 1.5])
	set(a, 'linewidth', 2);
	zeroAxes(gca, [], 10, 'arial')
	if exist('printDiagram') && printDiagram
	  fileName = ['dem' capName num2str(experimentNo)];
	  printPlot(fileName, '../tex/diagrams', '../html');
	end


	"""Shows that there is an optimum for the covariance function length scale.
	
	Description:
		DESC shows that by varying the length scale an artificial data
		set has different likelihoods, yet there is an optimum for which
		the likelihood is maximised.

	"""
		% COPYRIGHT : Neil D. Lawrence, 2006, 2008
	
	
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	fillColor = [0.7 0.7 0.7];
	markerSize = 20;
	markerWidth = 2;
	markerType = 'k.';
	lineWidth = 2;
	
	
	x = linspace(-1, 1, 6)';
	trueKern = kernCreate(x, {'rbf', 'white'});
	kern.comp{2}.variance = 0.001;
	K = kernCompute(trueKern, x);
	y = gsamp(zeros(1, 6), K, 1)';
	
	xtest = linspace(-1.5, 1.5, 200)';
	kern = trueKern;
	
	lengthScale = [0.05 0.1 0.25 0.5 1 2 4 8 16];
	counter = 0;
	
	figure(1)
	p = plot(x, y, markerType);
	set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 18)
	set(gca, 'ylim', [-2 1])
	set(gca, 'xlim', [-1.5 1.5])
	
	zeroAxes(gca);
	fileName = ['demOptimiseGp' num2str(counter)];
	if exist('printDiagram') && printDiagram
	  printPlot(fileName, '../tex/diagrams', '../html');
	end
	clf
	
	void = semilogx(NaN, NaN, 'k.-');
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 18)
	set(gca, 'ylim', [-12 -4])
	set(gca, 'xlim', [0.025 32]) 
	grid on
	ylabel('log-likelihood')
	xlabel('length scale')
	fileName = ['demOptimiseGp' num2str(counter) '0'];
	if exist('printDiagram') && printDiagram
	  printPlot(fileName, '../tex/diagrams', '../html');
	end
	
	clf
	
	for i = 1:length(lengthScale)
	  kern.comp{1}.inverseWidth = 1/(lengthScale(i)*lengthScale(i));
	  K = kernCompute(kern, x);
	  [invK, U] = pdinv(K);
	  logDetK = logdet(K, U);
	  ll(i) = -0.5*(logDetK + y'*invK*y + size(y, 1)*log(2*pi));
	  llLogDet(i) = -.5*(logDetK+size(y, 1)*log(2*pi));
	  llFit(i) = -.5*y'*invK*y;
	  Kx = kernCompute(kern, x, xtest);
	  ypredMean = Kx'*invK*y;
	  ypredVar = kernDiagCompute(kern, xtest) - sum((Kx'*invK).*Kx', 2);
	
	  counter = counter + 1;
	  figure(counter)
	  clf
	  fill([xtest; xtest(end:-1:1)], ...
	       [ypredMean; ypredMean(end:-1:1)] ...
	         + 2*[ypredVar; -ypredVar], ...
	       fillColor,'EdgeColor',fillColor)
	  hold on;
	  t = plot(xtest, ypredMean, 'k-');
	  
	  p = plot(x, y, markerType);
	  set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	  set(t, 'linewidth', lineWidth);
	  set(gca, 'fontname', 'times')
	  set(gca, 'fontsize', 18)
	  set(gca, 'ylim', [-2 1])
	  
	  zeroAxes(gca);
	  fileName = ['demOptimiseGp' num2str(counter)];
	  if exist('printDiagram') && printDiagram
	    printPlot(fileName, '../tex/diagrams', '../html');
	  end
	
	  counter = counter + 1;
	  figure(counter)
	  t = semilogx(lengthScale(1:i), ll(1:i), 'k.-');
	  hold on
	  t = [t; semilogx(lengthScale(1:i), llLogDet(1:i), 'k.:')];
	  t = [t; semilogx(lengthScale(1:i), llFit(1:i), 'k.--')];
	  set(t, 'markersize', markerSize, 'lineWidth', markerWidth);
	  set(gca, 'fontname', 'times')
	  set(gca, 'fontsize', 18)
	  set(gca, 'ylim', [-15 5])
	  set(gca, 'xlim', [0.025 32]) 
	  grid on
	  ylabel('log-likelihood')
	  xlabel('length scale')
	  fileName = ['demOptimiseGp' num2str(counter)];
	  if exist('printDiagram') && printDiagram
	    printPlot(fileName, '../tex/diagrams', '../html');
	  end
	end
	  
	
	

def model = gpSubspaceCreate(q,d,X,y,options,dim):

	"""
	
	Description:
		


	Copyright (c) 2008 Carl Henrik Ek
	
	"""
		  
	model = gpCreate(q,d,X,y,options);
	model.dim = dim;
	model.type = 'gpSubspace';
	
	return;

	"""Shows that there is an optimum for the covariance function length scale.
	
	Description:
		DESC shows that by varying the length scale an artificial data
		set has different likelihoods, yet there is an optimum for which
		the likelihood is maximised.

	"""
		% COPYRIGHT : Neil D. Lawrence, 2006, 2008
	
	
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	fillColor = [0.7 0.7 0.7];
	markerSize = 40;
	markerWidth = 4;
	markerType = 'k.';
	lineWidth = 4;
	
	
	x = linspace(-1, 1, 6)';
	trueKern = kernCreate(x, {'rbf', 'white'});
	kern.comp{2}.variance = 0.001;
	K = kernCompute(trueKern, x);
	y = gsamp(zeros(1, 6), K, 1)';
	
	xtest = linspace(-2, 2, 200)';
	kern = trueKern;
	
	lengthScale = [0.05 0.2 0.5 1 2 10];
	counter = 0;
	
	figure(1)
	p = plot(x, y, markerType);
	set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 24)
	set(gca, 'ylim', [-2 1])
	set(gca, 'xlim', [-2 2])
	
	zeroAxes(gca);
	
	void = semilogx(NaN, NaN, 'k.-');
	set(gca, 'fontname', 'times')
	set(gca, 'fontsize', 24)
	set(gca, 'ylim', [-12 -4])
	set(gca, 'xlim', [0.025 32]) 
	grid on
	ylabel('log-likelihood')
	xlabel('length scale')
	
	clf
	
	for i = 1:length(lengthScale)
	  kern.comp{1}.inverseWidth = 1/(lengthScale(i)*lengthScale(i));
	  K = kernCompute(kern, x);
	  [invK, U] = pdinv(K);
	  logDetK = logdet(K, U);
	  ll(i) = -0.5*(logDetK + y'*invK*y + size(y, 1)*log(2*pi));
	  llLogDet(i) = -.5*(logDetK+size(y, 1)*log(2*pi));
	  llFit(i) = -.5*y'*invK*y;
	  Kx = kernCompute(kern, x, xtest);
	  ypredMean = Kx'*invK*y;
	  ypredVar = kernDiagCompute(kern, xtest) - sum((Kx'*invK).*Kx', 2);
	
	  counter = counter + 1;
	  figure(counter)
	  clf
	  fill([xtest; xtest(end:-1:1)], ...
	       [ypredMean; ypredMean(end:-1:1)] ...
	         + 2*[ypredVar; -ypredVar], ...
	       fillColor,'EdgeColor',fillColor)
	  hold on;
	  t = plot(xtest, ypredMean, 'k-');
	  
	  p = plot(x, y, markerType);
	  set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	  set(t, 'linewidth', lineWidth);
	  set(gca, 'fontname', 'times')
	  set(gca, 'fontsize', 40)
	  set(gca, 'ylim', [-2 1])
	  
	  zeroAxes(gca, 0.025, 40);
	  fileName = ['demOptimiseGpTutorial' num2str(counter)];
	  if exist('printDiagram') && printDiagram
	    printPlot(fileName, '../tex/diagrams');
	  end
	  if i == length(lengthScale)
	    counter = counter + 1;
	    figure(counter)
	    t = semilogx(lengthScale(1:i), -ll(1:i), 'k.-');
	    hold on
	    t = [t; semilogx(lengthScale(1:i), -llLogDet(1:i), 'k.:')];
	    t = [t; semilogx(lengthScale(1:i), -llFit(1:i), 'k.--')];
	    set(t, 'markersize', markerSize/2, 'lineWidth', markerWidth/2);
	    set(gca, 'fontname', 'times')
	    set(gca, 'fontsize', 18)
	    set(gca, 'ylim', [0 15])
	    set(gca, 'xlim', [0.025 32]) 
	    grid on
	    ylabel('negative log-likelihood')
	    xlabel('length scale')
	    fileName = ['demOptimiseGpTutorial' num2str(counter)];
	    if exist('printDiagram') && printDiagram
	      printPlot(fileName, '../tex/diagrams');
	    end
	  end
	end
	  
	
	

def f = gpPosteriorSample(kernType, numSamps, params, lims, seed, bw):

	"""Create a plot of samples from a posterior covariance.
	
	Description:
	
	gpPosteriorSample(kernType, numSamps, params, lims) creates a plot
	 of samples from a kernel with the given parameters and variance.
	 Arguments:
	  kernType - the type of kernel to sample from.
	  numSamps - the number of samples to take.
	  params - parameter vector for the kernel.
	  lims - limits of the x axis.


	Copyright (c) 2008 Neil D. Lawrence
	
	"""
		
	global printDiagram
	  
	  
	if nargin < 6
	  bw = false;
	  if nargin < 5
	    seed = [];
	    if nargin < 4
	      lims = [-3 3];
	      if nargin < 3
	        params = [];
	        if nargin < 2
	          numSamps = 10;
	        end
	      end
	    end
	  end
	end
	if ~isempty(seed)
	  randn('seed', seed);
	  rand('seed', seed);
	end
	t_star = linspace(lims(1), lims(2), 200)';
	
	kern = kernCreate(t_star, kernType);
	for i=1:length(kern.comp)
	  kern.comp{i}.transforms = [];
	end
	
	if ~isempty(params)
	  feval = str2func([kern.type 'KernExpandParam']);
	  kern = feval(kern, params);
	end
	feval = str2func([kern.type 'KernExtractParam']);
	[params, names] = feval(kern);
	paramStr = [];
	for i = 1:length(names)
	  Name = names{i};
	  Name(1) = upper(Name(1));
	  ind = find(Name==' ');
	  Name(ind+1) = upper(Name(ind+1));
	  Name(ind) = '';
	  paramStr = [paramStr Name num2str(params(i))];
	  
	end
	paramStr(find(paramStr==46)) = 'p';
	infoStr = ['Samples' num2str(numSamps) 'Seed' num2str(randn('seed'))];
	
	% Covariance of the prior.
	K_starStar = kernCompute(kern, t_star, t_star);
	
	% Generate "training data" from a sine wave.
	t = rand(5, 1)*(lims(2)-lims(1))*0.6+lims(1)*0.6;
	f = sin(t);
	
	% Compute kernel for training data.
	K_starf = kernCompute(kern, t_star, t);
	K_ff = kernCompute(kern, t);
	
	% Mean and covariance of posterior.
	fbar = K_starf*pdinv(K_ff)*f;
	Sigma = K_starStar - K_starf*pdinv(K_ff)*K_starf';
	
	% Sample from the posterior.
	fsamp = real(gsamp(fbar, Sigma, numSamps));
	
	% Plot and save
	figure
	linHand = plot(t_star, fsamp);
	hold on
	linHandPlot = plot(t, f, 'r.')
	set(linHandPlot, 'markersize', 30)
	
	zeroAxes(gca, 0.025, 18, 'times');
	set(linHand, 'linewidth', 1)
	app = '';
	if bw
	  set(linHand, 'color', [0 0 0])
	  set(linHandPlot, 'color', [0 0 0])
	  app = 'bw';
	end
	
	if iscell(kernType)
	  KernType = [];
	  for i = length(kernType):-1:1
	    KernType = [kernType{i} KernType];
	    KernType(1) = upper(KernType(1));
	  end
	else
	  KernType(1) = upper(kernType(1));
	end
	
	if exist('printDiagram', 'var') & printDiagram
	  printPlot(['gpPosteriorSample' KernType infoStr paramStr app], ...
	            '../tex/diagrams', '../html')
	end
	
	


	"""Plot the true poses for the silhouette data.
	
	Description:
	
	"""
		dataSetName = 'silhouette';
	
	% load data
	[X, y, XTest, yTest] = mapLoadData(dataSetName);
	
	counter = 0;
	if printDiagram
	  fileBase = ['dem' capName 'GpTrue'];
	  for i = ind
	    counter = counter + 1;
	    figure
	    handle = xyzankurVisualise(yTest(i,:), 1);
	    printPlot([fileBase '_' num2str(counter)], '../tex/diagrams', '../html') 
	  end
	end


	"""Demonstrate Gaussian processes for regression on stick man data.
	
	Description:
	
	demStickGp1 runs a simple regression on the stick man data.
		

	See also
	GPCREATE, DEMINTERPOLATION


	Copyright (c) 2008 Neil D. Lawrence
	
	"""
		
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	dataSetName = 'stick';
	experimentNo = 1;
	
	% load data
	[y, lbls] = lvmLoadData(dataSetName);
	
	% load connectivity matrix
	[void, connect] = mocapLoadTextData([datasetsDirectory 'run1']);
	
	% Data is downsampled from 120 frames per second.
	fps = 120/4;
	t = (0:size(y, 1)-1)'/fps;
	% Predict Left Ankle X, Y and Z
	outputIndex = [7 8 9];
	indices = [1:10 20:55];
	testIndex = 1:size(y, 1);
	testIndex(indices) = [];
	yTrain = y(indices, outputIndex);
	tTrain = t(indices, 1);
	% Set up the model
	options = gpOptions('ftc');
	
	% Scale outputs to variance 1.
	options.scale2var1 = true;
	
	% Use the full Gaussian process model.
	q = size(tTrain, 2);
	d = size(yTrain, 2);
	model = gpCreate(q, d, tTrain, yTrain, options);
	
	display = 1;
	iters = 1000;
	
	model = gpOptimise(model, display, iters);
	modelDisplay(model)
	
	% Save results
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	save(['dem' capName 'Gp' num2str(experimentNo) '.mat'], 'model');
	
	
	% Plot results
	fillColor = [0.7 0.7 0.7];
	tTest = linspace(0, 2, 200)';
	[mu, varSigma] = gpPosteriorMeanVar(model, tTest);
	
	for i = 1:length(outputIndex);
	  figure
	  fill([tTest; tTest(end:-1:1)], ...
	       [mu(:, i); mu(end:-1:1, i)] ...
	       + 2*[sqrt(varSigma(:, i)); -sqrt(varSigma(end:-1:1, i))], ...
	       fillColor,'EdgeColor',fillColor)
	  hold on;
	  plot(t(testIndex), y(testIndex, outputIndex(i)), 'ko');
	  b=plot(tTrain, yTrain(:, i), 'k.');
	  set(b, 'Markersize', 15);
	  a = plot(tTest, mu(:, i), 'k-');
	  set(gca, 'xlim', [0 2])
	  set(a, 'linewidth', 2);
	  
	  zeroAxes(gca, [], 10, 'arial')
	  if exist('printDiagram') && printDiagram
	    fileName = ['dem' capName 'Gp' num2str(experimentNo) 'Out' num2str(i)];
	    printPlot(fileName, '../tex/diagrams', '../html');
	  end
	end
def model = gpUpdateAD(model, X):

	"""Update the representations of A and D associated with the model.
	
	Description:
	
	model = gpUpdateAD(model, X) updates the representations of A and
	 D in the model when called by gpUpdateKernels.
	 Returns:
	  model - the model with the A and D representations updated.
	 Arguments:
	  model - the model for which the representations are being updated.
	  X - the X values for which the representations are being computed.
		

	See also
	GPUPDATEKERNELS, GPEXPANDPARAM


	Copyright (c) 2009 Neil D. Lawrence 2006
	
	"""
		
	if nargin < 2
	  X = model.X;
	end
	
	switch model.approx
	 case 'ftc'
	  % Compute the inner product values.
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    for i = 1:model.d
	      model.innerProducts(1, i) = model.m(:, i)'*model.invK_uu...
	          *model.m(:, i);
	    end
	  else
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      model.innerProducts(1, i) = model.m(ind, i)'*model.invK_uu{i}...
	          *model.m(ind, i);
	    end
	  end
	  
	 case {'dtc', 'dtcvar'}
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    % Compute A = invBetaK_uu + K_uf*K_uf'
	    K_uf2 = model.K_uf*model.K_uf';
	    model.A = (1/model.beta)*model.K_uu+ K_uf2;
	    % This can become unstable when K_uf2 is low rank.
	    [model.Ainv, U] = pdinv(model.A);
	    model.logdetA = logdet(model.A, U);
	 
	    % compute inner products
	    for i = 1:model.d
	      E = model.K_uf*model.m(:, i);    
	      model.innerProducts(1, i) = ...
	          model.beta*(model.m(:, i)'*model.m(:, i) ...
	                      - E'*model.Ainv*E);
	    end
	    if strcmp(model.approx, 'dtcvar')
	      model.diagD = model.beta*(model.diagK ...
	        - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)');
	    end
	  else
	    if ~model.isMissingData
	      K_uf2 = model.K_uf*model.K_uf';
	    end
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      % Compute A = invBetaK_uu + K_uf*K_uf'
	      if model.isMissingData
	        K_uf2 = model.K_uf(:, ind)*model.K_uf(:, ind)';
	      end
	      model.A{i} = (1/model.beta)*model.K_uu+ K_uf2;
	      % This can become unstable when K_uf2 is low rank.
	      [model.Ainv{i}, U] = pdinv(model.A{i});
	      model.logdetA(i) = logdet(model.A{i}, U);
	      % compute inner products
	      E = model.K_uf(:, ind)*model.m(ind, i);    
	      model.innerProducts(1, i) = ...
	          model.beta*(model.m(ind, i)'*model.m(ind, i) ...
	                      - E'*model.Ainv{i}*E);
	    end
	    if strcmp(model.approx, 'dtcvar')
	      error('Non spherical implementation for dtcvar not yet done.')
	    end
	  end
	  
	 case 'fitc'
	  model.L = jitChol(model.K_uu)';
	
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    model.diagD = 1 + model.beta*model.diagK ...
	        - model.beta*sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
	    % model.diagDdivBeta = 1/model.beta + model.diagK - sum(model.K_uf.*(model.invK_uu*model.K_uf), 1)';
	    model.Dinv = sparseDiag(1./model.diagD);
	    % model.betaDinv = sparseDiag(1./model.diagDdivBeta);
	    K_ufDinvK_uf = model.K_uf*model.Dinv*model.K_uf';
	    model.A = 1/model.beta*model.K_uu + K_ufDinvK_uf;    
	    % model.betaA = model.K_uu 
	    % This can become unstable when K_ufDinvK_uf is low rank.
	    [model.Ainv, U] = pdinv(model.A);
	    model.logDetA = logdet(model.A, U);
	    % model.detDiff = model.logDetA - model.logDetK_uu;
	    model.detDiff = - log(model.beta)*model.k + log(det(eye(model.k) + model.beta*K_ufDinvK_uf*model.invK_uu));
	    % compute inner products
	    for i = 1:model.d
	      Dinvm = model.Dinv*model.m(:, i);
	      K_ufDinvm = model.K_uf*Dinvm;
	      model.innerProducts(1, i) = model.beta*(Dinvm'*model.m(:, i) ...
	          - K_ufDinvm'*model.Ainv*K_ufDinvm);
	    end
	  
	    % Computations from Ed's implementation.
	    model.V = model.L\model.K_uf;
	    model.V = model.V./repmat(sqrt(model.diagD)', model.k, 1);
	    model.Am = 1/model.beta*eye(model.k)+model.V*model.V';
	    model.Lm = jitChol(model.Am)';
	    model.invLmV = model.Lm\model.V;
	    model.scaledM = model.m./repmat(sqrt(model.diagD), 1, model.d);
	    model.bet = model.invLmV*model.scaledM;
	  else
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      model.diagD{i} = 1 + model.beta*model.diagK(ind) ...
	          - model.beta*sum(model.K_uf(:, ind).*(model.invK_uu*model.K_uf(:, ind)), 1)';
	      model.Dinv{i} = sparseDiag(1./model.diagD{i});
	      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
	          *model.K_uf(:, ind)';
	      model.A{i} = 1/model.beta*model.K_uu + K_ufDinvK_uf;
	      % This can become unstable when K_ufDinvK_uf is low rank.
	      [model.Ainv{i}, U] = pdinv(model.A{i});
	      model.logDetA(i) = logdet(model.A{i}, U);
	      % model.detDiff = model.logDetA - model.logDetK_uu;
	      model.detDiff(i) = - log(model.beta)*model.k + log(det(eye(model.k) + model.beta*K_ufDinvK_uf*model.invK_uu));
	    
	      % compute inner products
	      Dinvm = model.Dinv{i}*model.m(ind, i);
	      K_ufDinvm = model.K_uf(:, ind)*Dinvm;
	      model.innerProducts(1, i) = model.beta*(Dinvm'*model.m(ind, i) - K_ufDinvm'*model.Ainv{i}*K_ufDinvm);
	    
	      % Computations from Ed's implementation.
	      model.V{i} = model.L\model.K_uf(:, ind);
	      model.V{i} = model.V{i}./repmat(sqrt(model.diagD{i})', model.k, 1);
	      model.Am{i} = 1/model.beta*eye(model.k)+model.V{i}*model.V{i}';
	      model.Lm{i} = jitChol(model.Am{i})';
	      model.invLmV{i} = model.Lm{i}\model.V{i};
	      model.scaledM{i} = model.m(ind, i)./sqrt(model.diagD{i});
	      model.bet{i} = model.invLmV{i}*model.scaledM{i};
	
	    end
	      
	  end
	  
	 case 'pitc'
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    model.A = 1/model.beta*model.K_uu;
	    K_ufDinvm = zeros(model.k, model.d);
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      model.D{i} = eye(length(ind)) + model.beta*model.K{i} - ...
	          model.beta*model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
	      [model.Dinv{i}, U] = pdinv(model.D{i});
	      model.logDetD(i) = logdet(model.D{i}, U);
	      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}...
	          *model.K_uf(:, ind)';
	      model.A = model.A + K_ufDinvK_uf;
	      Dinvm{i} = model.Dinv{i}*model.m(ind, :);
	      K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
	    end
	    % This can become unstable when K_ufDinvK_uf is low rank.
	    [model.Ainv, U] = pdinv(model.A);
	    model.logDetA = logdet(model.A, U);
	    % compute inner products
	    for i = 1:model.d
	      model.innerProducts(1, i) = - model.beta*K_ufDinvm(:, i)'*model.Ainv*K_ufDinvm(:, ...
	                                                        i);
	    end
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      for j = 1:model.d
	        model.innerProducts(1, j) = model.innerProducts(1, j) ...
	            + model.beta*Dinvm{i}(:, j)'*model.m(ind, j);
	      end
	    end
	  else
	    for j = 1:model.d
	      model.A{j} = 1/model.beta*model.K_uu;
	      K_ufDinvm = zeros(model.k, model.d);
	      for i = 1:length(model.blockEnd)
	        ind = gpDataIndices(model, j, i);
	        model.D{i, j} = eye(length(ind)) + model.beta*model.K{i, j} - ...
	            model.beta*model.K_uf(:, ind)'*model.invK_uu*model.K_uf(:, ind);
	        [model.Dinv{i, j}, U] = pdinv(model.D{i, j});
	        model.logDetD(i, j) = logdet(model.D{i, j}, U);
	        K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i, j}...
	            *model.K_uf(:, ind)';
	        model.A{j} = model.A{j} + K_ufDinvK_uf;
	        Dinvm{i}(ind, j) = model.Dinv{i, j}*model.m(ind, j);
	        K_ufDinvm(:, j) = K_ufDinvm(:, j) + model.K_uf(:, ind)*Dinvm{i}(ind, ...
	                                                          j);
	      end
	      % This can become unstable when K_ufDinvK_uf is low rank.
	      [model.Ainv{j}, U] = pdinv(model.A{j});
	      model.logDetA(j) = logdet(model.A{j}, U);
	    end
	    
	    % compute inner products
	    for j = 1:model.d
	      model.innerProducts(1, j) = - model.beta*K_ufDinvm(:, j)'*model.Ainv{j}*K_ufDinvm(:, ...
	                                                        j);
	    end
	    for i = 1:length(model.blockEnd)
	      for j = 1:model.d
	        ind = gpDataIndices(model, j, i);
	        model.innerProducts(1, j) = model.innerProducts(1, j) ...
	            + model.beta*Dinvm{i}(ind, j)'*model.m(ind, j);
	      end
	    end
	    
	  end
	 otherwise
	  error('Unknown approximating criterion.')
	end
	

def [f, g] = gpObjectiveGradient(params, model):

	"""Wrapper function for GP objective and gradient.
	
	Description:
	
	[f, g] = gpObjectiveGradient(params, model) returns the negative
	 log likelihood of a Gaussian process model given the model
	 structure and a vector of parameters. This allows the use of
	 NETLAB minimisation functions to find the model parameters.
	 Returns:
	  f - the negative log likelihood of the GP model.
	  g - the gradient of the negative log likelihood of the GP model
	   with respect to the parameters.
	 Arguments:
	  params - the parameters of the model for which the objective will
	   be evaluated.
	  model - the model structure for which the objective will be
	   evaluated.
		

	See also
	MINIMIZE, GPCREATE, GPGRADIENT, GPLOGLIKELIHOOD, GPOPTIMISE


	Copyright (c) 2005, 2006 Neil D. Lawrence
	
	"""
		
	% Check how the optimiser has given the parameters
	if size(params, 1) > size(params, 2)
	  % As a column vector ... transpose everything.
	  transpose = true;
	  model = gpExpandParam(model, params');
	else
	  transpose = false;
	  model = gpExpandParam(model, params);
	end
	
	f = - gpLogLikelihood(model);
	if nargout > 1
	  g = - gpLogLikeGradients(model);
	end
	if transpose
	  g = g';
	end

	"""Do a simple 1-D regression after Snelson & Ghahramani's example.
	
	Description:
	
	"""
		
	% Fix seeds
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	dataSetName = 'spgp1d';
	experimentNo = 4;
	
	% load data
	[X, y] = mapLoadData(dataSetName);
	
	% Set up the model
	options = gpOptions('ftc');
	options.optimiser = 'conjgrad';
	
	% Make use of correct kernel.
	options.kern = kernCreate(X, {'rbf', 'white'});
	options.kern.comp{1}.inverseWidth = 20;
	options.kern.comp{2}.variance = 0.01;
	
	% Use the full Gaussian process model.
	q = size(X, 2);
	d = size(y, 2);
	model = gpCreate(q, d, X, y, options);
	
	% This updates the kernels.
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);
	
	% Save results
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	save(['dem' capName num2str(experimentNo) '.mat'], 'model');
	
	
	demSpgp1dPlot
def g = gpGradient(params, model):

	"""Gradient wrapper for a GP model.
	
	Description:
	
	g = gpGradient(params, model) wraps the log likelihood gradient
	 function to return the gradient of the negative of the log
	 likelihood. This can then be used in, for example, NETLAB,
	 minimisation tools.
	 Returns:
	  g - the returned gradient of the negative log likelihood for the
	   given parameters.
	 Arguments:
	  params - the parameters of the model.
	  model - the model for which gradients will be computed.
		

	See also
	SCG, CONJGRAD, GPCREATE, GPGRADIENT, GPLOGLIKEGRADIENT, GPOPTIMISE


	Copyright (c) 2005, 2006 Neil D. Lawrence
	
	"""
		
	model = gpExpandParam(model, params);
	g = - gpLogLikeGradients(model);

def model = gpExpandParam(model, params):

	"""Expand a parameter vector into a GP model.
	
	Description:
	
	model = gpExpandParam(model, params) takes the given vector of
	 parameters and places them in the model structure, it then updates
	 any stored representations that are dependent on those parameters,
	 for example kernel matrices etc..
	 Returns:
	  model - a returned model structure containing the updated
	   parameters.
	 Arguments:
	  model - the model structure for which parameters are to be
	   updated.
	  params - a vector of parameters for placing in the model
	   structure.
		

	See also
	GPCREATE, GPEXTRACTPARAM, MODELEXTRACTPARAM, GPUPDATEKERNELS


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	if strcmp(model.approx, 'ftc') | model.fixInducing
	  endVal = 0;
	else
	  startVal = 1;
	  endVal = model.k*model.q;
	  model.X_u = reshape(params(startVal:endVal), model.k, model.q);
	end
	startVal = endVal +1;
	endVal = endVal + model.kern.nParams;
	model.kern = kernExpandParam(model.kern, params(startVal:endVal));
	
	% Check if there is a mean function.
	if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
	  startVal = endVal + 1;
	  endVal = endVal + model.meanFunction.numParams;
	  model.meanFunction = modelExpandParam(model.meanFunction, ...
	                                        params(startVal:endVal));
	end
	
	% Check if the output scales are being learnt.
	if model.learnScales
	  startVal = endVal + 1;
	  endVal = endVal + model.d;
	  fhandle = str2func([model.scaleTransform 'Transform']);
	  model.scale = fhandle(params(startVal:endVal), 'atox');
	  model.m = gpComputeM(model);
	end
	
	% Check if beta is being optimised.
	if model.optimiseBeta
	  startVal = endVal + 1;
	  endVal = endVal + prod(size(model.beta));
	  fhandle = str2func([model.betaTransform 'Transform']);
	  model.beta = fhandle(params(startVal:endVal), 'atox');
	end
	
	% Record the total number of parameters.
	model.nParams = endVal;
	
	% Update the kernel representations.
	switch model.approx
	 case 'ftc'
	  model = gpUpdateKernels(model, model.X, model.X_u);
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  model = gpUpdateKernels(model, model.X, model.X_u);
	 otherwise
	  error('Unknown approximation type.')
	end
	
	% Update the vector 'alpha' for computing posterior mean.
	if isfield(model, 'alpha')
	  model = gpComputeAlpha(model);
	end

def gpWriteToFile(fileName, model):

	"""Write a file to be read by the C++ implementation.
	
	Description:
	
	gpWriteToFile(fileName, model) writes a GP model to a file
	 produced by the C++ GP implementation.
	 Arguments:
	  fileName - the file name written by the C++ software.
	  model - a MATLAB GP model structure containing the model for the
	   file.
		

	See also
	GPWRITETOFID, GPCREATE


	Copyright (c) 2008 Neil D. Lawrence
	
	"""
		
	FID = fopen(fileName, 'w');
	if FID==-1
	  error(['Cannot open file ' fileName])
	end
	gpWriteToFID(FID, model);
	fclose(FID);
def y = gpOut(model, x):

	"""Evaluate the output of an Gaussian process model.
	
	Description:
	
	y = gpOut(model, x) evaluates the output of a given Gaussian
	 process model.
	 Returns:
	  y - the output of the GP model. The function checks if there is a
	   noise model associated with the GP, if there is, it is used,
	   otherwise the mean of gpPosteriorMeanVar is returned.
	 Arguments:
	  model - the model for which the output is being evaluated.
	  x - the input position for which the output is required.
		

	See also
	GPCREATE, GPPOSTERIORMEANVAR


	Copyright (c) 2006 Neil D. Lawrence and Carl Ek
	
	"""
		
	if nargin < 2
	  % This implies evaluate for the training data.
	  mu = model.mu;
	  varsigma = model.varSigma;
	else
	  if isfield(model, 'noise')
	    [mu, varsigma] = gpPosteriorMeanVar(model, x);
	    y = noiseOut(model.noise, mu, varsigma);
	  else
	    y = gpPosteriorMeanVar(model, x);
	  end
	end


	"""Shows the average of the poses.
	
	Description:
	
	"""
		% FORMAT
	% DESC show the average of the training poses from Agarwal and Triggs.
	%
	% SEEALSO : demSilhouetteGp1, demSilhouetteGp2
	% 
	% COPYRIGHT : Neil D. Lawrence, 2008
	
	
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	dataSetName = 'silhouette';
	
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	
	% load data
	[X, y, XTest, yTest] = mapLoadData(dataSetName);
	
	yPred = repmat(mean(y), size(yTest, 1), 1);
	xyzankurAnimCompare(yPred, yTest);
	
	yDiff = (yPred - yTest);
	rmsError = sqrt(sum(sum(yDiff.*yDiff))/prod(size(yDiff)));
	
	counter = 0;
	if printDiagram
	  fileBase = ['dem' capName 'Average'];
	  figure
	  handle = xyzankurVisualise(yPred(1, :), 1);
	  printPlot([fileBase], '../tex/diagrams', '../html') 
	end


	"""Do a simple 1-D regression after Snelson & Ghahramani's example.
	
	Description:
	
	"""
		
	% Fix seeds
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	dataSetName = 'spgp1d';
	experimentNo = 1;
	
	% load data
	[X, y] = mapLoadData(dataSetName);
	
	% Set up model
	options = gpOptions('dtc');
	options.numActive = 9;
	options.optimiser = 'conjgrad';
	
	% use the deterministic training conditional.
	q = size(X, 2);
	d = size(y, 2);
	
	model = gpCreate(q, d, X, y, options);
	model.X_u = randn(9, 1)*0.25 - 0.75;
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);
	
	% Optimise the model.
	iters = 1000;
	display = 1;
	
	model = gpOptimise(model, display, iters);
	
	% Save the results.
	capName = dataSetName;
	capName(1) = upper(capName(1));
	save(['dem' capName num2str(experimentNo) '.mat'], 'model');
	
	
	demSpgp1dPlot


	"""Demonstrate Gaussian processes for regression.
	
	Description:
	
	demRegression runs a simple one-D Gaussian process displaying
	 errorbars.
		

	See also
	GPCREATE, DEMINTERPOLATION


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	noiseLevel = 0.2;
	noiseVar = noiseLevel*noiseLevel;
	% Create data set
	x = linspace(-1, 1, 9)';
	trueKern = kernCreate(x, 'rbf');
	K = kernCompute(trueKern, x) + eye(size(x, 1))*noiseVar;
	% Sample some true function values.
	yTrue = gsamp(zeros(size(x))', K, 1)';
	
	fillColor =[0.7 0.7 0.7];
	markerSize = 20;
	markerWidth = 6;
	markerType = 'k.';
	lineWidth = 2;
	% Create a test set
	indTrain{1} = [1 9]';
	indTrain{2} = [1 5 9]';
	indTrain{3} = [1 3 5 7 9]';
	indTrain{4} = [1 2 3 4 5 6 7 8 9]';
	figNo = 1;
	for i = 0:length(indTrain)
	  if i > 0
	    yTrain = yTrue(indTrain{i});
	    xTrain = x(indTrain{i});
	    kern = kernCreate(x, 'rbf');
	    % Change inverse variance (1/(lengthScale^2)))
	    kern.inverseWidth = 5;
	    
	    xTest = linspace(-2, 2, 200)';
	    
	    Kx = kernCompute(kern, xTest, xTrain);
	    Ktrain = kernCompute(kern, xTrain, xTrain);
	    
	    yPred = Kx*pdinv(Ktrain + eye(size(Ktrain))*noiseVar)*yTrain;
	    yVar = kernDiagCompute(kern, xTest) - sum(Kx*pdinv(Ktrain+ eye(size(Ktrain))*noiseVar).*Kx, 2);
	    ySd = sqrt(yVar);
	    figure(figNo)
	    clf
	    fill([xTest; xTest(end:-1:1)], ...
	         [yPred; yPred(end:-1:1)] ...
	         + 2*[ySd; -ySd], ...
	         fillColor,'EdgeColor',fillColor)
	    hold on;
	    h = plot(xTest, yPred, 'k-');
	    set(h, 'linewidth', lineWidth)
	    p = plot(xTrain, yTrain, markerType);
	    set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
	    set(gca, 'xtick', [-2 -1 0 1 2]);
	    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
	    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
	    zeroAxes(gca);
	    if exist('printDiagram') && printDiagram
	      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
	    end
	    figNo = figNo + 1;
	  else
	    xTest = linspace(-2, 2, 200)';
	    kern = kernCreate(xTest, 'rbf');
	    % Change inverse variance (1/(lengthScale^2)))
	    kern.inverseWidth = 5;
	    
	    figure(figNo)
	    p = [];
	    yPred = zeros(size(xTest));
	    ySd = sqrt(kernDiagCompute(kern, xTest));
	    fill([xTest; xTest(end:-1:1)], ...
	         [yPred; yPred(end:-1:1)] ...
	         + 2*[ySd; -ySd], ...
	         fillColor,'EdgeColor',fillColor)
	    hold on;
	    h = plot(xTest, yPred, 'k-');
	    
	    set(h, 'linewidth', lineWidth)
	    set(gca, 'xtick', [-2 -1 0 1 2]);
	    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
	    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
	    zeroAxes(gca);
	    if exist('printDiagram') && printDiagram
	      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
	    end
	    figNo = figNo + 1;
	    
	  end
	  if i < length(indTrain)
	    figure(figNo)
	    clf
	    fill([xTest; xTest(end:-1:1)], ...
	         [yPred; yPred(end:-1:1)] ...
	         + 2*[ySd; -ySd], ...
	         fillColor,'EdgeColor',fillColor)
	    hold on;
	    h = plot(xTest, yPred, 'k-');
	    set(h, 'linewidth', lineWidth)
	    p = [p plot(x(indTrain{i+1}), yTrue(indTrain{i+1}), markerType)];
	    set(p, 'markersize', markerSize, 'linewidth', markerWidth);
	    set(gca, 'xtick', [-2 -1 0 1 2]);
	    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
	    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
	    zeroAxes(gca);
	    if exist('printDiagram') && printDiagram
	      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
	    end
	    figNo = figNo + 1;
	  end
	end
def ll = gpLogLikelihood(model):

	"""Compute the log likelihood of a GP.
	
	Description:
	
	ll = gpLogLikelihood(model) computes the log likelihood of a data
	 set given a GP model.
	 Returns:
	  ll - the log likelihood of the data in the GP model.
	 Arguments:
	  model - the GP model for which log likelihood is to be computed.
		

	See also
	GPCREATE, GPLOGLIKEGRADIENTS, MODELLOGLIKELIHOOD


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	switch model.approx
	 case 'ftc'
	  % No approximation, just do a full computation on K.
	  ll = 0;
	  for i = 1:size(model.m, 2)
	    if ~isfield(model, 'isSpherical') | model.isSpherical
	      ll = ll -.5*model.logDetK_uu- .5*model.m(:, i)'*model.invK_uu*model.m(:, i);
	    else
	      if model.isMissingData
	        m = model.m(model.indexPresent{i}, i);
	      else
	        m = model.m(:, i);
	      end
	      ll = ll - .5*model.logDetK_uu(i) - .5*m'*model.invK_uu{i}*m;
	    end
	  end
	 case {'dtc', 'dtcvar'}
	  % Deterministic training conditional
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    E = model.K_uf*model.m;
	    EET = E*E';
	    if length(model.beta)==1
	      ll =  -0.5*(model.d*(-(model.N-model.k)*log(model.beta) ...
	                           - model.logDetK_uu +model.logdetA) ...
	                  - (sum(sum(model.Ainv.*EET)) ...
	                     -sum(sum(model.m.*model.m)))*model.beta);
	      if strcmp(model.approx, 'dtcvar')
	        ll = ll - model.d*0.5*sum(model.diagD);
	      end
	    else
	      error('Not implemented variable length beta yet.');
	    end
	  else
	    ll = 0;
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      e = model.K_uf(:, ind)*model.m(ind, i);
	      if length(model.beta)==1
	        ll = ll - 0.5*((-(model.N-model.k)*log(model.beta) ...
	                        - model.logDetK_uu +model.logdetA(i)) ...
	                       - (e'*model.Ainv{i}*e ...
	                          -model.m(ind, i)'*model.m(ind, i))* ...
	                       model.beta);
	        if(isnan(ll))
	          error('Log likelihood is NaN')
	        end
	        if strcmp(model.approx, 'dtcvar')
	          error('Not implemented dtcvar for non-spherical yet.');
	        end
	      else
	        error('Not implemented variable length beta yet.');
	      end
	    end
	  end
	 case 'fitc'
	  % Fully independent training conditional.
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    if length(model.beta)==1
	      if false
	        % This is the original objective
	        Dinvm = model.Dinv*model.m;
	        K_ufDinvm = model.K_uf*Dinvm;
	        ll = -0.5*(model.d*(sum(log(model.diagD))...
	                            -(model.N-model.k)*log(model.beta) ...
	                            + model.detDiff)...
	                   + (sum(sum(Dinvm.*model.m))...
	                      - sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm)))*model.beta);
	        
	        ll = ll - 0.5*model.N*model.d*log(2*pi);
	      else
	        % This is objective to match Ed Snelson's code
	        ll =  - model.d*(sum(log(diag(model.Lm))) + 0.5*(-(model.N - model.k)*log(model.beta)+(model.N*log(2*pi)+sum(log(model.diagD)))));
	        for i = 1:model.d
	          ll = ll - 0.5*model.beta*(model.scaledM(:, i)'*model.scaledM(:, i) ...
	                                    - model.bet(:, i)'*model.bet(:, i));
	        end
	      end
	    else
	      error('Variable length Beta not implemented yet.')
	    end
	  else
	    if length(model.beta)==1
	      if false
	        ll = 0;
	        for i = 1:model.d
	          ind = gpDataIndices(model, i);
	          Dinvm = model.Dinv{i}*model.m(ind, i);
	          K_ufDinvm = model.K_uf(:, ind)*Dinvm;
	          ll = ll -0.5*(sum(log(model.diagD{i})) ...
	                        - (length(ind) - model.k)*log(model.beta) ...
	                        + model.detDiff(i) ...
	                        + (sum(sum(Dinvm.*model.m(ind, i))) ...
	                           - sum(sum((model.Ainv{i}*K_ufDinvm).* ...
	                                     K_ufDinvm)))*model.beta ...
	                        +length(ind)*log(2*pi));
	        end
	      else
	        % This is objective to match Ed Snelson's code
	        ll = 0;
	        for i = 1:model.d
	          ind = gpDataIndices(model, i);
	          ll =  ll - (sum(log(diag(model.Lm{i}))) ...
	                      + 0.5*(-(length(ind) - model.k)*log(model.beta) ...
	                             +(length(ind)*log(2*pi)+sum(log(model.diagD{i})))));
	          ll = ll - 0.5*model.beta*(model.scaledM{i}'*model.scaledM{i} ...
	                                    - model.bet{i}'*model.bet{i});
	        end
	      end
	    else
	      error('Variable length Beta not implemented yet.')
	    end
	  end
	 case 'pitc'
	  % Partially independent training conditional.
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    if length(model.beta)==1
	      ll = model.d*(model.logDetA-model.logDetK_uu +model.k*log(model.beta));
	      % Loop through the blocks computing each part to be added.
	      K_ufDinvm = zeros(model.k, model.d);
	      for i = 1:length(model.blockEnd)
	        ind = gpBlockIndices(model, i);
	        Dinvm{i} = model.Dinv{i}*model.m(ind, :);
	        K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
	      end
	      ll = ll - model.beta*sum(sum((model.Ainv*K_ufDinvm).*K_ufDinvm));
	      
	      for i = 1:length(model.blockEnd)
	        ind = gpBlockIndices(model, i); 
	        ll = ll + model.d*(model.logDetD(i) ...
	                           - length(ind)*log(model.beta))...
	             + model.beta*sum(sum(Dinvm{i}.*model.m(ind, :)));
	      end
	      ll = -0.5*ll;
	      ll = ll - 0.5*model.N*model.d*log(2*pi);
	    else
	      error('Variable Length Beta not implemented yet.')
	    end
	  else
	    if length(model.beta)==1
	      
	      ll = 0;
	      for j = 1:model.d
	        ll = ll + model.logDetA(j)-model.logDetK_uu + model.k*log(model.beta);
	        % Loop through the blocks computing each part to be added.
	        K_ufDinvm = zeros(model.k, 1);
	        for i = 1:length(model.blockEnd)
	          ind = gpDataIndices(model, j, i);
	          Dinvm{i, j} = model.Dinv{i, j}*model.m(ind, j);
	          K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i, j};
	        end
	        ll = ll - model.beta*sum(sum((model.Ainv{i}*K_ufDinvm).*K_ufDinvm));
	        
	        for i = 1:length(model.blockEnd)
	          ind = gpDataIndices(model, j, i);
	          ll = ll + model.logDetD(i, j) ...
	               - length(ind)*log(model.beta) ...
	               + model.beta*sum(sum(Dinvm{i, j}.*model.m(ind, j)));
	          ll = ll + length(ind)*log(2*pi);
	        end
	      end
	      ll = -0.5*ll;
	    else
	      error('Variable Length Beta not implemented yet.');
	    end
	  end
	end
	if model.learnScales
	  ll = ll - sum(log(model.scale));
	end
	ll = ll - model.d*model.N/2*log(2*pi);
def f = gpSample(kernType, numSamps, params, lims, seed, bw):

	"""Create a plot of samples from a GP.
	
	Description:
	
	gpSample(kernType, numSamps, params, lims) creates a plot of
	 samples from a kernel with the given parameters and variance.
	 Arguments:
	  kernType - the type of kernel to sample from.
	  numSamps - the number of samples to take.
	  params - parameter vector for the kernel.
	  lims - limits of the x axis.


	Copyright (c) 2008 Neil D. Lawrence
	
	"""
		  
	if nargin < 6
	  bw = false;
	  if nargin < 5
	    seed = [];
	    if nargin < 4
	      lims = [-3 3];
	      if nargin < 3
	        params = [];
	        if nargin < 2
	          numSamps = 10;
	        end
	      end
	    end
	  end
	end
	if ~isempty(seed)
	  randn('seed', seed);
	  rand('seed', seed);
	end
	t_star = linspace(lims(1), lims(2), 200)';
	
	kern = kernCreate(t_star, kernType);
	if ~isempty(params)
	  feval = str2func([kernType 'KernExpandParam']);
	  kern = feval(kern, params);
	end
	feval = str2func([kernType 'KernExtractParam']);
	[params, names] = feval(kern);
	paramStr = [];
	for i = 1:length(names)
	  Name = names{i};
	  Name(1) = upper(Name(1));
	  ind = find(Name==' ');
	  Name(ind+1) = upper(Name(ind+1));
	  Name(ind) = '';
	  paramStr = [paramStr Name num2str(params(i))];
	  
	end
	infoStr = ['Samples' num2str(numSamps) 'Seed' num2str(randn('seed'))];
	
	% Covariance of the prior.
	K_starStar = kernCompute(kern, t_star, t_star);
	
	% Sample from the prior.
	fsamp = real(gsamp(zeros(size(t_star)), K_starStar, numSamps));
	
	% Plot and save.
	clf
	linHand = plot(t_star, fsamp);
	zeroAxes(gca, 0.025, 18, 'times');
	set(linHand, 'linewidth', 2)
	if bw
	  set(linHand, 'color', [0 0 0])
	end
	KernType = kernType;
	KernType(1) = upper(kernType(1));
	fileName = ['gpSample' KernType infoStr paramStr];
	printPlot(fileName, '../tex/diagrams', '../html');
	


	"""Model silhouette data with independent MLP GPs.
	
	Description:
	
	"""
		% FORMAT
	% DESC runs a simple regression on the Agawal and Triggs data.
	%
	% SEEALSO : gpCreate, demInterpolation
	% 
	% COPYRIGHT : Neil D. Lawrence, 2008
	
	
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	dataSetName = 'silhouette';
	experimentNo = 2;
	
	% load data
	[X, y, XTest, yTest] = mapLoadData(dataSetName);
	
	% Set up the model
	options = gpOptions('ftc');
	options.kern{1} = 'mlp';
	
	% Scale outputs to variance 1.
	options.scale2var1 = true;
	
	% Use the full Gaussian process model.
	q = size(X, 2);
	d = size(y, 2);
	model = gpCreate(q, d, X, y, options);
	
	display = 1;
	iters = 1000;
	
	model = gpOptimise(model, display, iters);
	modelDisplay(model)
	
	% Save results
	capName = dataSetName;
	capName(1) = upper(capName(1));
	fileBaseName = ['dem' capName 'Gp' num2str(experimentNo)];
	save([fileBaseName '.mat'], 'model');
	demSilhouettePlot

def ind = gpDataIndices(model, dimNo, blockNo):

	"""Return indices of present data.
	
	Description:
	
	ind = gpDataIndices(model, dimNo) returns the indices of data
	 which is not missing for a given dimension in the GP-LVM.
	 Returns:
	  ind - indices of training data along that dimension which isn't
	   missing.
	 Arguments:
	  model - the model for which the indices are being returned.
	  dimNo - the dimension for which the presence of missing data is
	   being looked at.
		DESC returns the indices of data which is not missing for a given
		dimension in the GP-LVM and a block number in the PITC approximation.
		ARG model : the model for which the indices are being returned.
		ARG dimNo : the dimension for which the presence of missing data
		is being looked at.
		ARG blockNo : the block number in the PITC approximation for
		which the indices are required.
		RETURN ind : indices of training data along that dimension which
		isn't missing.
		
		

	See also
	GPCREATE


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	
	if nargin > 2
	  if ~strcmp(model.approx, 'pitc')
	    error('Block number only relevant for pitc approximation.');
	  else
	    if blockNo == 1
	      startVal = 1;
	    else
	      startVal = model.blockEnd(blockNo-1)+1;
	    end
	    endVal = model.blockEnd(blockNo);
	    if model.isMissingData
	      st = min(find(model.indexPresent{dimNo}>=startVal));
	      fi = max(find(model.indexPresent{dimNo}<=endVal));
	      ind = model.indexPresent{dimNo}(st:fi)';
	    else
	      ind = startVal:endVal;
	    end
	  end
	else
	  if strcmp(model.approx, 'pitc')
	    error('Must give block number with PITC approximation');
	  else
	    if model.isMissingData
	      ind = model.indexPresent{dimNo}';
	    else
	      ind = 1:model.N;
	    end
	  end
	end

def [gParam, gX_u, gX, g_beta] = gpLogLikeGradients(model, X, M, X_u):

	"""Compute the gradients for the parameters and X.
	
	Description:
	
	gParam = gpLogLikeGradients(ARG) computes the gradients of the
	 Gaussian process log likelihood with respect to the parameters of
	 the model.
	 Returns:
	  gParam - the gradient of the log likelihood with respect to the
	   model parameters.
	 Arguments:
	  ARG - the model structure for which gradients are computed.
		DESC computes the gradients of the Gaussian process log
		likelihood with respect to the parameters of the model and with
		respect to any inducing variables.
		ARG : model : the model structure for which gradients are computed.
		RETURN gParam : the gradient of the log likelihood with respect to
		the model parameters.
		RETURN gX_u : the gradient of the log likelihood with respect to
		the inducing variables. If inducing variables aren't being used
		this returns zero.
		
		DESC computes the gradients of the Gaussian process log
		likelihood with respect to the parameters of the model, with
		respect to any inducing variables and with respect to input data
		locations. This is used for computing gradients in the GP-LVM.
		ARG : model : the model structure for which gradients are computed.
		RETURN gParam : the gradient of the log likelihood with respect to
		the model parameters (including any gradients with respect to beta).
		RETURN gX_u : the gradient of the log likelihood with respect to
		the inducing variables. If inducing variables aren't being used
		this returns zero.
		RETURN gX : the gradient of the log likelihood with respect to
		the input data locations.
		
		DESC computes the gradients of the Gaussian process log
		likelihood with respect to the parameters of the model, with
		respect to any inducing variables and with respect to input data
		locations. This is used for computing gradients in the GP-LVM.
		ARG : model : the model structure for which gradients are computed.
		RETURN gParam : the gradient of the log likelihood with respect to
		the model parameters.
		RETURN gX_u : the gradient of the log likelihood with respect to
		the inducing variables. If inducing variables aren't being used
		this returns zero.
		RETURN gX : the gradient of the log likelihood with respect to
		the input data locations.
		RETURN gbeta : the gradient of the log likelihood with respect to beta.
		
		DESC computes the gradients of the Gaussian process log
		likelihood with respect to the model parameters (and optionally,
		as above with respect to inducing variables and input data) given
		the target data, input data and inducing variable
		locations.
		ARG : model : the model structure for which gradients are computed.
		ARG : X : the input data locations for which gradients are computed.
		ARG : M : the scaled and bias removed target data for which the
		gradients are computed.
		ARG : X_U : the inducing variable locations for which gradients are computed.
		RETURN gParam : the gradient of the log likelihood with respect to
		the model parameters.
		
		
		

	See also
	GPLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS, FGPLVMLOGLIKEGRADIENTS


	Copyright (c) 2005, 2006, 2007, 2009 Neil D. Lawrence
	
	
	With modifications by Carl Henrik Ek 2008
	
	"""
		
	if nargin < 4
	  if isfield(model, 'X_u')
	    X_u = model.X_u;
	  else
	    X_u = [];
	  end
	  if nargin < 3
	    M = model.m;
	  end
	  if nargin < 2
	    X = model.X;
	  end
	end
	
	gX_u = [];
	gX = [];
	
	g_scaleBias = gpScaleBiasGradient(model);
	if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
	  g_meanFunc = gpMeanFunctionGradient(model);
	else
	  g_meanFunc = [];
	end
	
	switch model.approx
	 case 'ftc'
	  % Full training conditional.
	  if nargout > 2
	    %%% Prepare to Compute Gradients with respect to X %%%
	    gKX = kernGradX(model.kern, X, X);
	    gKX = gKX*2;
	    dgKX = kernDiagGradX(model.kern, X);
	    for i = 1:model.N
	      gKX(i, :, i) = dgKX(i, :);
	    end
	    gX = zeros(model.N, model.q);
	  end
	  
	  %%% Gradients of Kernel Parameters %%%
	  g_param = zeros(1, model.kern.nParams);
	  if isfield(model, 'beta')
	    g_beta = 0;
	  else
	    g_beta = [];
	  end
	  for k = 1:model.d
	    gK = localCovarianceGradients(model, M(:, k), k);
	    if nargout > 2
	      %%% Compute Gradients with respect to X %%%
	      ind = gpDataIndices(model, k);
	      counter = 0;
	      for i = ind
	        counter = counter + 1;
	        for j = 1:model.q
	          gX(i, j) = gX(i, j) + gKX(ind, j, i)'*gK(:, counter);
	        end
	      end
	    end
	    %%% Compute Gradients of Kernel Parameters %%%
	    if model.isMissingData
	      g_param = g_param ...
	               + kernGradient(model.kern, ...
	                              X(model.indexPresent{k}, :), ...
	                              gK);
	    else
	      g_param = g_param + kernGradient(model.kern, X, gK);
	    end
	    if isfield(model, 'beta') && model.optimiseBeta
	      if size(model.beta, 1) == 1
	        g_beta = g_beta + sum(diag(gK));
	      elseif size(model.beta, 2)==1 ...
	            & size(model.beta, 1)==model.N
	        g_beta = g_beta + diag(gK);
	      elseif size(model.beta, 2) == model.d ...
	            & size(model.beta, 1) == model.N
	        g_beta(:, k) = diag(gK);
	      else
	        error('Unusual dimensions for model.beta.');
	      end
	    end
	  end
	   
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  % Sparse approximations.
	  [gK_u, gK_uf, gK_star, g_beta] = gpCovGrads(model, M);
	  
	  %%% Compute Gradients of Kernel Parameters %%%
	  gParam_u = kernGradient(model.kern, X_u, gK_u);
	  gParam_uf = kernGradient(model.kern, X_u, X, gK_uf);
	
	  g_param = gParam_u + gParam_uf;
	  
	  %%% Compute Gradients with respect to X_u %%%
	  gKX = kernGradX(model.kern, X_u, X_u);
	  
	  % The 2 accounts for the fact that covGrad is symmetric
	  gKX = gKX*2;
	  dgKX = kernDiagGradX(model.kern, X_u);
	  for i = 1:model.k
	    gKX(i, :, i) = dgKX(i, :);
	  end
	  
	  if ~model.fixInducing | nargout > 1
	    % Allocate space for gX_u
	    gX_u = zeros(model.k, model.q);
	    % Compute portion associated with gK_u
	    for i = 1:model.k
	      for j = 1:model.q
	        gX_u(i, j) = gKX(:, j, i)'*gK_u(:, i);
	      end
	    end
	    
	    % Compute portion associated with gK_uf
	    gKX_uf = kernGradX(model.kern, X_u, X);
	    for i = 1:model.k
	      for j = 1:model.q
	        gX_u(i, j) = gX_u(i, j) + gKX_uf(:, j, i)'*gK_uf(i, :)';
	      end
	    end
	  end
	  if nargout > 2
	    %%% Compute gradients with respect to X %%%
	    
	    % Allocate space for gX
	    gX = zeros(model.N, model.q);
	    
	    % this needs to be recomputed so that it is wrt X not X_u
	    gKX_uf = kernGradX(model.kern, X, X_u);
	    
	    for i = 1:model.N
	      for j = 1:model.q
	        gX(i, j) = gKX_uf(:, j, i)'*gK_uf(:, i);
	      end
	    end    
	  end
	 otherwise
	  error('Unknown model approximation.')
	end
	
	
	switch model.approx
	 case 'ftc'
	  % Full training conditional. Nothing required here.
	 case 'dtc'
	  % Deterministic training conditional.  
	 case {'fitc', 'dtcvar'}
	  % Fully independent training conditional.
	  % Variational sparse approximation.
	  
	  if nargout > 2
	    % deal with diagonal term's effect on X gradients..
	    gKXdiag = kernDiagGradX(model.kern, X);
	    for i = 1:model.N
	      gX(i, :) = gX(i, :) + gKXdiag(i, :)*gK_star(i);
	    end
	  end
	  
	  % deal with diagonal term's affect on kernel parameters.
	  g_param = g_param + kernDiagGradient(model.kern, X, gK_star);
	
	
	 case 'pitc'
	  % Partially independent training conditional.
	  
	  if nargout > 2
	    % deal with block diagonal term's effect on X gradients.
	    startVal = 1;
	    for i = 1:length(model.blockEnd)
	      endVal = model.blockEnd(i);
	      ind = startVal:endVal;
	      gKXblock = kernGradX(model.kern, X(ind, :), X(ind, :));
	      
	      % The 2 accounts for the fact that covGrad is symmetric
	      gKXblock = gKXblock*2;
	      
	      % fix diagonal
	      dgKXblock = kernDiagGradX(model.kern, X(ind, :));
	      for j = 1:length(ind)
	        gKXblock(j, :, j) = dgKXblock(j, :);
	      end
	      
	      for j = ind
	        for k = 1:model.q
	          subInd = j - startVal + 1;
	          gX(j, k) = gX(j, k) + gKXblock(:, k, subInd)'*gK_star{i}(:, subInd);
	        end
	      end
	      startVal = endVal + 1;
	    end
	  end
	  % deal with block diagonal's effect on kernel parameters.
	  for i = 1:length(model.blockEnd);
	    ind = gpBlockIndices(model, i);
	    g_param = g_param ...
	              + kernGradient(model.kern, X(ind, :), gK_star{i});
	  end
	  
	 otherwise
	  error('Unrecognised model approximation');
	end
	
	if nargout < 4
	  if (~isfield(model, 'optimiseBeta') && ~strcmp(model.approx, 'ftc')) ...
	      | model.optimiseBeta
	    % append beta gradient to end of parameters
	    gParam = [g_param(:)' g_meanFunc g_scaleBias g_beta];
	  else
	    gParam = [g_param(:)' g_meanFunc g_scaleBias];
	  end
	else
	    gParam = [g_param(:)' g_meanFunc g_scaleBias];
	end
	
	% if there is only one output argument, pack gX_u and gParam into it.
	if nargout == 1;
	  gParam = [gX_u(:)' gParam];
	end
	
def gK = localCovarianceGradients(model, y, dimension):
	
	% LOCALCOVARIANCEGRADIENTS
	
	if ~isfield(model, 'isSpherical') || model.isSpherical
	  invKy = model.invK_uu*y;
	  gK = -model.invK_uu + invKy*invKy';
	else
	  if model.isMissingData
	    m = y(model.indexPresent{dimension});
	  else
	    m = y;
	  end
	  invKy = model.invK_uu{dimension}*m;
	  gK = -model.invK_uu{dimension} + invKy*invKy';
	end
	gK = gK*.5;
	    

def [gK_uu, gK_uf, g_Lambda, gBeta] = gpCovGrads(model, M):

	"""Sparse objective function gradients wrt Covariance functions for inducing variables.
	
	Description:
	
	[gK_uu, gK_uf, gLambda, gBeta] = gpCovGrads(model, M) gives the
	 gradients of the log likelihood with respect to the components of
	 the sparse covariance (or the full covariance for the ftc case).
	 Returns:
	  gK_uu - the gradient of the likelihood with respect to the
	   elements of K_uu (or in the case of the 'ftc' criterion the
	   gradients with respect to the kernel).
	  gK_uf - the gradient of the likelihood with respect to the
	   elements of K_uf.
	  gLambda - the gradient of the likelihood with respect to the
	   diagonal term in the fitc approximation and the blocks of the pitc
	   approximation.
	  gBeta - the gradient with respect to the beta term in the
	   covariance structure.
	 Arguments:
	  model - the model for which the gradients are to be computed.
	  M - The training data for which the computation is to be made
		

	See also
	GPCREATE, GPLOGLIKEGRADIENT


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	switch model.approx
	 case {'dtc', 'dtcvar'}
	  % Deterministic training conditional.
	  if strcmp(model.approx, 'dtcvar')
	    dtcvar = true;
	  else
	    dtcvar = false;
	  end
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    E = model.K_uf*M;
	    EET = E*E';
	    AinvEET = model.Ainv*EET;
	    AinvEETAinv = AinvEET*model.Ainv;
	    gK_uu = 0.5*(model.d*(model.invK_uu-(1/model.beta)*model.Ainv) ...
	                 - AinvEETAinv);
	    if dtcvar
	      K_uuInvK_uf = model.invK_uu*model.K_uf;
	      gK_uu = gK_uu - 0.5*model.d*model.beta...
	              *K_uuInvK_uf*K_uuInvK_uf';
	    end
	    AinvK_uf = model.Ainv*model.K_uf;
	    gK_uf = -model.d*AinvK_uf-model.beta*(AinvEET*AinvK_uf-(model.Ainv*E*M'));
	    if dtcvar
	      gK_uf = gK_uf + model.d*model.beta*K_uuInvK_uf;
	    end
	    gBeta = 0.5*(model.d*((model.N-model.k)/model.beta ...
	                              +sum(sum(model.Ainv.*model.K_uu))/(model.beta*model.beta))...
	                     +sum(sum(AinvEETAinv.*model.K_uu))/model.beta ...
	                     +(trace(AinvEET)-sum(sum(M.*M))));
	    if dtcvar
	      gBeta = gBeta -0.5*model.d*sum(model.diagD)/model.beta;
	    end
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	    if dtcvar
	      g_Lambda = repmat(-0.5*model.beta*model.d, 1, model.N);
	    else
	      g_Lambda = [];
	    end
	  else
	    gK_uu = zeros(model.k, model.k);
	    gK_uf = zeros(model.k, model.N);
	    gBeta = 0;
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      e = model.K_uf(:, ind)*M(ind, i);
	      Ainve = model.Ainv{i}*e;
	      AinveeT = Ainve*e';      
	      AinveeTAinv = Ainve*Ainve';
	      gK_uu = gK_uu+0.5*((model.invK_uu-(1/model.beta)*model.Ainv{i}) ...
	                         - AinveeTAinv);
	      
	      AinvK_uf = model.Ainv{i}*model.K_uf(:, ind);
	      gK_uf(:, ind) = gK_uf(:, ind) - AinvK_uf...
	          -model.beta*(AinveeT*AinvK_uf-(Ainve*M(ind, i)'));
	      
	      gBeta = gBeta ...
	          + 0.5*(((model.N-model.k)/model.beta ...
	                  +sum(sum(model.Ainv{i}.*model.K_uu))/(model.beta*model.beta))...
	                 +sum(sum(AinveeTAinv.*model.K_uu))/model.beta ...
	                 +(trace(AinveeT)-sum(sum(M(ind, i).*M(ind, i)))));
	    end
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	    g_Lambda = [];
	  end
	 case 'fitc'
	  % Fully independent training conditonal.
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    E = model.K_uf*model.Dinv*M;
	    EET = E*E';
	    AinvE = model.Ainv*E;
	    %AinvEET = model.Ainv*EET;
	    diagK_fuAinvEMT = sum(model.K_uf.*(model.Ainv*E*M'), 1)';
	    AinvEETAinv = AinvE*AinvE';
	    diagK_ufdAinvplusAinvEETAinvK_fu = ...
	        sum(model.K_uf.*((model.d*model.Ainv+model.beta*AinvEETAinv)*model.K_uf), 1)';
	    invK_uuK_uf = model.invK_uu*model.K_uf;
	    if true
	      invK_uuK_ufDinv = invK_uuK_uf*model.Dinv;
	    else
	      invK_uuK_ufDinv = model.L'\model.V;
	    end
	    diagMMT = sum(M.*M, 2);
	    diagQ = -model.d*model.diagD + model.beta*diagMMT ...
	            + diagK_ufdAinvplusAinvEETAinvK_fu...
	            -2*model.beta*diagK_fuAinvEMT;
	    gK_uu = 0.5*(model.d*(model.invK_uu ...
	                 -model.Ainv/model.beta) - AinvEETAinv ...
	                 + model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*invK_uuK_ufDinv');
	    gK_uf = -model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*model.Dinv ...      
	            -model.d*model.Ainv*model.K_uf*model.Dinv ...
	            -model.beta*AinvEETAinv*model.K_uf*model.Dinv ...
	            +model.beta*model.Ainv*E*M'*model.Dinv;
	    g_Lambda = (0.5*diagQ*model.beta)./(model.diagD.*model.diagD);
	    gBeta = -sum(g_Lambda)/(model.beta*model.beta);
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	  else
	    gK_uu = zeros(model.k, model.k);
	    gK_uf = zeros(model.k, model.N);
	    g_Lambda = zeros(model.N, 1);
	    gBeta = 0;
	    for i = 1:model.d
	      ind = gpDataIndices(model, i);
	      K_ufDinvK_uf = model.K_uf(:, ind)*model.Dinv{i}*model.K_uf(:, ...
	                                                        ind)';
	      e = model.K_uf(:, ind)*model.Dinv{i}*M(ind, i);
	      Ainve = model.Ainv{i}*e;
	      AinveeTAinv = Ainve*Ainve';
	      diagK_fuAinveyT = sum(model.K_uf(:, ind).*(Ainve*M(ind,i)'), 1)';
	      diagK_ufdAinvplusAinveeTAinvK_fu = ...
	          sum(model.K_uf(:, ind).*((model.Ainv{i}+model.beta*AinveeTAinv)*model.K_uf(:, ...
	                                                        ind)), 1)';
	      invK_uuK_uf = model.invK_uu*model.K_uf(:, ind);
	      invK_uuK_ufDinv = invK_uuK_uf*model.Dinv{i};
	      diagyyT = M(ind, i).*M(ind, i);
	      diagQ = -model.diagD{i} + model.beta*diagyyT ...
	              + diagK_ufdAinvplusAinveeTAinvK_fu...
	              -2*model.beta*diagK_fuAinveyT;
	      gK_uu = gK_uu ...
	              +0.5*(model.invK_uu ...
	                    - model.Ainv{i}/model.beta - AinveeTAinv ...
	                    + model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*invK_uuK_ufDinv');
	      gK_uf(:, ind) = gK_uf(:, ind) ...
	              -model.beta*invK_uuK_ufDinv*sparseDiag(diagQ)*model.Dinv{i} ...      
	              -model.Ainv{i}*model.K_uf(:, ind)*model.Dinv{i} ...
	              -model.beta*AinveeTAinv*model.K_uf(:, ind)*model.Dinv{i} ...
	              +model.beta*Ainve*M(ind, i)'*model.Dinv{i};
	      g_Lambda(ind) = g_Lambda(ind) ...
	          + 0.5*model.beta*diagQ./(model.diagD{i}.*model.diagD{i});
	    end
	    gBeta = gBeta - sum(g_Lambda)/(model.beta*model.beta);
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	  end
	 
	 case 'pitc' 
	  % Partially independent training conditional.
	  if ~isfield(model, 'isSpherical') | model.isSpherical
	    E = zeros(model.k, model.d);
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      E = E + model.K_uf(:, ind)*model.Dinv{i}*M(ind, :);
	    end
	    AinvE = model.Ainv*E;
	    AinvEET = AinvE*E';
	    AinvEETAinv = AinvEET*model.Ainv;
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      K_fuAinvEMT = model.beta*model.K_uf(:, ind)'*AinvE*M(ind, :)';
	      blockQ{i} = -model.d*model.D{i} + model.beta*M(ind, :)*M(ind, :)' ...
	          + model.K_uf(:, ind)'*(model.d*model.Ainv + model.beta*AinvEETAinv)*model.K_uf(:, ind)...
	          -K_fuAinvEMT - K_fuAinvEMT';
	    end
	    gK_uu = model.d*model.invK_uu ...
	            - model.d*model.Ainv/model.beta - AinvEETAinv;
	    gBeta = 0;
	    gK_ufBase = -(model.d*model.Ainv + model.beta*AinvEETAinv)*model.K_uf ...
	        + model.beta*AinvE*M';
	    
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      invK_uuK_ufDinv = model.invK_uu*model.K_uf(:, ind)*model.Dinv{i};
	      gK_uu = gK_uu + model.beta*invK_uuK_ufDinv*blockQ{i}*invK_uuK_ufDinv';
	      
	      gK_uf(:, ind) = (gK_ufBase(:, ind) ...
	                       -model.beta*invK_uuK_ufDinv*blockQ{i})*model.Dinv{i};
	      
	      
	      g_Lambda{i} = 0.5*model.Dinv{i}*blockQ{i}*model.Dinv{i}*model.beta;
	      gBeta = gBeta - sum(diag((g_Lambda{i})))/(model.beta*model.beta);
	    end
	    gK_uu = gK_uu*0.5;
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	  else
	    gK_uu = zeros(model.k, model.k);
	    gK_uf = zeros(model.k, model.N);
	    for i = 1:length(model.blockEnd)
	      if i == 1
	        indLen = model.blockEnd(1);
	      else
	        indLen = model.blockEnd(i) - model.blockEnd(i-1);
	      end
	      g_Lambda{i} = zeros(indLen, indLen);
	    end
	    gBeta = 0;
	    for j = 1:model.d
	      e = zeros(model.k, 1);
	      for i = 1:length(model.blockEnd)
	        ind = gpDataIndices(model, j, i);
	        e = e + model.K_uf(:, ind)*model.Dinv{i, j}*M(ind, j);
	      end
	      Ainve = model.Ainv{j}*e;
	      AinveeT = Ainve*e';
	      AinveeTAinv = AinveeT*model.Ainv{j};
	      for i = 1:length(model.blockEnd)
	        ind = gpDataIndices(model, j, i);
	        K_fuAinveyT = model.beta*model.K_uf(:, ind)'*Ainve*M(ind, j)';
	        blockQ{i} = -model.D{i, j} + model.beta*M(ind, j)*M(ind, j)' ...
	            + model.K_uf(:, ind)'*(model.Ainv{j} + model.beta*AinveeTAinv)*model.K_uf(:, ind)...
	            -K_fuAinveyT - K_fuAinveyT';
	      end
	      gK_uu = gK_uu + model.invK_uu ...
	            - model.Ainv{j}/model.beta - AinveeTAinv;
	    
	      for i = 1:length(model.blockEnd)
	        ind = gpDataIndices(model, j, i);
	        gK_ufBase = -(model.Ainv{i} + model.beta*AinveeTAinv)*model.K_uf(:, ind) ...
	            + model.beta*Ainve*M(ind, j)';
	        invK_uuK_ufDinv = model.invK_uu*model.K_uf(:, ind)*model.Dinv{i,j};
	        gK_uu = gK_uu + model.beta*invK_uuK_ufDinv*blockQ{i}*invK_uuK_ufDinv';
	      
	        gK_uf(:, ind) = gK_uf(:, ind) + (gK_ufBase ...
	                         -model.beta*invK_uuK_ufDinv*blockQ{i})*model.Dinv{i,j};
	      
	        if i == 1
	          localInd = ind;
	        else
	          localInd = ind - (model.blockEnd(i-1));
	        end
	        g_Lambda{i}(localInd, localInd) = g_Lambda{i}(localInd, localInd) ...
	            + 0.5*model.Dinv{i,j}*blockQ{i}*model.Dinv{i,j}*model.beta;
	      end
	    end
	    for i = 1:length(model.blockEnd)
	      gBeta = gBeta - sum(diag((g_Lambda{i})))/(model.beta*model.beta);
	    end
	    gK_uu = gK_uu*0.5;
	    fhandle = str2func([model.betaTransform 'Transform']);
	    gBeta = gBeta*fhandle(model.beta, 'gradfact');
	  end
	 otherwise
	  error('Unknown approximation type');
	end

	"""Do a simple 1-D regression after Snelson & Ghahramani's example.
	
	Description:
	
	"""
		
	% Fix seeds
	randn('seed', 1e5);
	rand('seed', 1e5);
	
	dataSetName = 'spgp1d';
	experimentNo = 3;
	
	% load data
	[X, y] = mapLoadData(dataSetName);
	
	% Set up model
	options = gpOptions('pitc');
	options.numActive = 9;
	options.optimiser = 'conjgrad';
	
	% use the deterministic training conditional.
	q = size(X, 2);
	d = size(y, 2);
	
	model = gpCreate(q, d, X, y, options);
	model.X_u = randn(9, 1)*0.25 - 0.75;
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);
	
	% Optimise the model.
	iters = 1000;
	display = 1;
	
	model = gpOptimise(model, display, iters);
	
	% Save the results.
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	save(['dem' capName num2str(experimentNo) '.mat'], 'model');
	
	demSpgp1dPlot
def [kern, noise, gpInfo] = gpDeconstruct(model):

	"""break GP in pieces for saving.
	
	Description:
	
	[kern, noise, gpInfo] = gpDeconstruct(model) takes an GP model
	 structure and breaks it into component parts for saving.
	 Returns:
	  kern - the kernel component of the GP model.
	  noise - the noise component of the GP model.
	  gpInfo - a structure containing the other information from the GP:
	   what the sparse approximation is, what the inducing variables are.
	 Arguments:
	  model - the model that needs to be saved.
		

	See also
	GPRECONSTRUCT


	Copyright (c) 2007, 2009 Neil D. Lawrence
	
	"""
		
	kern = model.kern;
	if isfield(model, 'noise')
	  noise = model.noise;
	else
	  noise = [];
	end
	gpInfo.learnScales = model.learnScales;
	gpInfo.approx = model.approx;
	switch model.approx
	 case 'ftc'
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  gpInfo.beta = model.beta;
	  gpInfo.betaTransform = model.betaTransform;
	  gpInfo.fixInducing = model.fixInducing;
	  if model.fixInducing
	    gpInfo.inducingIndices = model.inducingIndices;
	  else
	    gpInfo.X_u = model.X_u;
	  end
	end
	gpInfo.type = 'gp';
	gpInfo.scale = model.scale;
	gpInfo.bias = model.bias;
	gpInfo.d = model.d;
	gpInfo.q = model.q;
	gpInfo.k = model.k;
	gpInfo.N = model.N;
def model = gpCreate(q, d, X, y, options):

	"""Create a GP model with inducing varibles/pseudo-inputs.
	
	Description:
	
	model = gpCreate(q, d, X, y, options) creates a Gaussian process
	 model structure with default parameter settings as specified by
	 the options vector.
	 Returns:
	  model - model structure containing the Gaussian process.
	 Arguments:
	  q - input data dimension.
	  d - the number of processes (i.e. output data dimension).
	  X - the input data matrix.
	  y - the target (output) data.
	  options - options structure as defined by gpOptions.m.
		
		

	See also
	GPOPTIONS, MODELCREATE


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	
	With modifications by Cark Henrik Ek 2007
	
	"""
		
	if size(X, 2) ~= q
	  error(['Input matrix X does not have dimension ' num2str(q)]);
	end
	if size(y, 2) ~= d
	  error(['Input matrix y does not have dimension ' num2str(d)]);
	end
	
	if any(isnan(y)) & ~options.isMissingData
	  error('NaN values in y, but no missing data declared.')
	end
	if options.isMissingData & options.isSpherical
	  error('If there is missing data, spherical flag cannot be set.');
	end
	
	model.type = 'gp';
	model.approx = options.approx;
	  
	model.learnScales = options.learnScales;
	model.scaleTransform = optimiDefaultConstraint('positive');
	
	model.optimiseBeta = options.optimiseBeta;
	model.betaTransform =  optimiDefaultConstraint('positive');  
	
	model.X = X;
	model.y = y;
	
	model.q = size(X, 2);
	model.d = size(y, 2);
	model.N = size(y, 1);
	
	% Set up a mean function if one is given.
	if isfield(options, 'meanFunction') & ~isempty(options.meanFunction)
	  if isstruct(options.meanFunction)
	    model.meanFunction = options.meanFunction;
	  else
	    if ~isempty(options.meanFunction)
	      model.meanFunction = modelCreate(options.meanFunction, model.q, model.d, options.meanFunctionOptions);
	    end
	  end
	end
	
	
	model.optimiser = options.optimiser;
	
	model.isMissingData = options.isMissingData;
	if model.isMissingData
	  for i = 1:model.d
	    model.indexPresent{i} = find(~isnan(y(:, i)));
	  end
	end
	
	model.isSpherical = options.isSpherical;
	
	if ~model.isMissingData
	  model.bias = mean(y);
	  model.scale = ones(1, model.d);
	else
	  for i = 1:model.d
	    if isempty(model.indexPresent{i})
	      model.bias(i) = 0;
	      model.scale(i) = 1;
	    else
	      model.bias(i) = mean(model.y(model.indexPresent{i}, i));
	      model.scale(i) = 1;
	    end
	  end
	end
	if(isfield(options,'scale2var1'))
	  if(options.scale2var1)
	    model.scale = std(model.y);
	    model.scale(find(model.scale==0)) = 1;
	    if(model.learnScales)
	      warning('Both learn scales and scale2var1 set for GP');
	    end
	  end
	end
	
	
	model.m = gpComputeM(model);
	
	if isstruct(options.kern) 
	  model.kern = options.kern;
	else
	  model.kern = kernCreate(model.X, options.kern);
	end
	
	if isfield(options, 'noise')
	  if isstruct(options.noise)
	    model.noise = options.noise;
	  else
	    model.noise = noiseCreate(options.noise, y);
	  end
	
	  % Set up noise model gradient storage.
	  model.nu = zeros(size(y));
	  model.g = zeros(size(y));
	  model.gamma = zeros(size(y));
	  
	  % Initate noise model
	  model.noise = noiseCreate(noiseType, y); 
	  
	  % Set up storage for the expectations
	  model.expectations.f = model.y;
	  model.expectations.ff = ones(size(model.y));
	  model.expectations.fBar =ones(size(model.y));
	  model.expectations.fBarfBar = ones(numData, ...
	                                     numData, ...
	                                     size(model.y, 2));
	end
	
	
	
	switch options.approx
	 case 'ftc'
	  model.k = 0;
	  model.X_u = [];
	  if model.optimiseBeta
	    model.beta = options.beta
	    if isempty(options.beta)
	      error('options.beta cannot be empty if it is being optimised.');
	    end
	  end
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  % Sub-sample inducing variables.
	  model.k = options.numActive;
	  model.fixInducing = options.fixInducing;
	  if options.fixInducing
	    if length(options.fixIndices)~=options.numActive
	      error(['Length of indices for fixed inducing variables must ' ...
	             'match number of inducing variables']);
	    end
	    model.X_u = model.X(options.fixIndices, :);
	    model.inducingIndices = options.fixIndices;
	  else
	    ind = randperm(model.N);
	    ind = ind(1:model.k);
	    model.X_u = model.X(ind, :);
	  end
	  model.beta = options.beta;
	end
	if model.k>model.N
	  error('Number of active points cannot be greater than number of data.')
	end
	if strcmp(model.approx, 'pitc')
	  numBlocks = ceil(model.N/model.k);
	  numPerBlock = ceil(model.N/numBlocks);
	  startVal = 1;
	  endVal = model.k;
	  model.blockEnd = zeros(1, numBlocks);
	  for i = 1:numBlocks
	    model.blockEnd(i) = endVal;
	    endVal = numPerBlock + endVal;
	    if endVal>model.N
	      endVal = model.N;
	    end
	  end  
	end
	
	initParams = gpExtractParam(model);
	
	% This forces kernel computation.
	model = gpExpandParam(model, initParams);
	
	

def model = gpSubspaceExpandParam(model,params):

	"""
	
	Description:
		


	Copyright (c) 2008 Carl Henrik Ek
	
	"""
		  
	model = gpExpandParam(model,params);
	
	return;
def [params, names] = gpExtractParam(model):

	"""Extract a parameter vector from a GP model.
	
	Description:
	
	params = gpExtractParam(model) extracts the model parameters from
	 a structure containing the information about a Gaussian process.
	 Returns:
	  params - a vector of parameters from the model.
	 Arguments:
	  model - the model structure containing the information about the
	   model.
		DESC does the same as above, but also returns parameter names.
		ARG model : the model structure containing the information about
		the model.
		RETURN params : a vector of parameters from the model.
		RETURN names : cell array of parameter names.
		
		

	See also
	GPCREATE, GPEXPANDPARAM, MODELEXTRACTPARAM


	Copyright (c) 2005, 2006, 2007, 2009 Neil D. Lawrence
	
	"""
		
	if nargout > 1
	  returnNames = true;
	else
	  returnNames = false;
	end
	
	% Check if the output scales are being learnt.
	if model.learnScales
	  fhandle = str2func([model.scaleTransform 'Transform']);
	  scaleParams = fhandle(model.scale, 'xtoa');
	  if returnNames
	    for i = 1:length(scaleParams)
	      scaleParamNames{i} = ['Output Scale ' num2str(i)];
	    end
	  end
	else
	  scaleParams = [];
	  scaleParamNames = {};
	end
	
	% Check if there is a mean function.
	if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
	  if returnNames
	    [meanFuncParams, meanFuncParamNames] = modelExtractParam(model.meanFunction);
	    for i = 1:length(meanFuncParamNames)
	      meanFuncParamNames{i} = ['Mean Func, ' meanFuncParamNames{i}];
	    end
	  else
	    meanFuncParams = modelExtractParam(model.meanFunction);
	  end
	else
	  meanFuncParamNames = {};
	  meanFuncParams =[];
	end
	if returnNames
	  [kernParams, kernParamNames] = kernExtractParam(model.kern);
	  for i = 1:length(kernParamNames)
	    kernParamNames{i} = ['Kernel, ' kernParamNames{i}];
	  end
	else
	  kernParams = kernExtractParam(model.kern);
	end
	switch model.approx
	 case 'ftc'
	  params =  [kernParams meanFuncParams scaleParams];
	  if returnNames
	    names = {kernParamNames{:}, meanFuncParamNames{:}, scaleParamNames{:}};
	  end
	  if model.optimiseBeta
	    fhandle = str2func([model.betaTransform 'Transform']);
	    betaParam = fhandle(model.beta, 'xtoa');
	    params = [params betaParam(:)'];
	    if returnNames
	      for i = 1:length(betaParam)
	        betaParamNames{i} = ['Beta ' num2str(i)];
	      end
	      names = {names{:}, betaParamNames{:}};
	    end
	  end
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  paramPart = [kernParams meanFuncParams scaleParams];
	  if returnNames
	    names = {kernParamNames{:}, meanFuncParamNames{:}, ...
	             scaleParamNames{:}};
	  end
	  if model.optimiseBeta
	    fhandle = str2func([model.betaTransform 'Transform']);
	    betaParam = fhandle(model.beta, 'xtoa');
	    paramPart = [paramPart betaParam(:)'];
	    if returnNames
	      for i = 1:length(betaParam)
	        betaParamNames{i} = ['Beta ' num2str(i)];
	      end
	      names = {names{:}, betaParamNames{:}};
	    end
	  end
	  if model.fixInducing
	    params = paramPart;
	  else
	    params =  [model.X_u(:)' paramPart];
	    if returnNames 
	      for i = 1:size(model.X_u, 1)
	        for j = 1:size(model.X_u, 2)
	          X_uNames{i, j} = ['X_u(' num2str(i) ', ' num2str(j) ')'];
	        end
	      end
	      names = {X_uNames{:}, names{:}};
	    end
	  end
	end
def model = gpSubspaceOptimise(model,varargin):

	""
	
	Description:
		


	Copyright (c) 2008 Carl Henrik Ek
	
	"""
		
	model = gpOptimise(model,varargin{:});
	
	return;
def options = gpOptions(approx):

	"""Return default options for GP model.
	
	Description:
	
	options = gpOptions(approx) returns the default options in a
	 structure for a GP model.
	 Returns:
	  options - structure containing the default options for the given
	   approximation type.
	 Arguments:
	  approx - approximation type, either 'ftc' (no approximation),
	   'dtcvar' (variational sparse approximation) 'dtc' (deterministic
	   training conditional), 'fitc' (fully independent training
	   conditional) or 'pitc' (partially independent training
	   conditional).
		

	See also
	GPCREATE


	Copyright (c) 2005, 2006, 2007, 2009 Neil D. Lawrence
	
	"""
		
	if nargin < 1
	  options.approx = 'ftc';
	else
	  options.approx = approx;
	end
	
	% Select type of optimiser.
	options.optimiser = optimiDefaultOptimiser;
	
	% Set to true to learn output scales.
	options.learnScales = false;
	
	% Set to true to scale outputs to variance 1.
	options.scale2var1 = false;
	
	% Set to true to optimise beta.
	switch approx
	 case 'ftc'
	  options.optimiseBeta = false;
	 otherwise
	  options.optimiseBeta = true;
	end
	
	% Set to a given mean function to have a mean function.
	options.meanFunction = [];
	% Options structure for mean function options.
	options.meanFunctionOptions = [];
	
	% Set to 1 if output processes have a shared variance.
	options.isSpherical = 1;
	
	% Set to 1 if there is data missing in the target matrix.
	options.isMissingData = 0;
	
	switch options.approx
	 case 'ftc'
	  % bog-standard kernel.
	  options.kern = {'rbf', 'bias', 'white'};
	  options.numActive = 0;
	  options.beta = [];
	 case {'fitc', 'pitc', 'dtc', 'dtcvar'}
	  options.kern = {'rbf', 'bias', 'white'};
	  options.numActive = 100;
	  options.beta = 1e3;
	
	  % Option to fix the inducing variables to other latent points.
	  options.fixInducing = 0;
	  options.fixIndices = [];
	end
	

def model = gpReconstruct(kern, noise, gpInfo, X, y):

	"""Reconstruct an GP form component parts.
	
	Description:
	
	model = gpReconstruct(kern, noise, gpInfo, X, y) takes component
	 parts of an GP model and reconstructs the GP model. The component
	 parts are normally retrieved from a saved file.
	 Returns:
	  model - an GP model structure that combines the component parts.
	 Arguments:
	  kern - a kernel structure for the GP.
	  noise - a noise structure for the GP (currently ignored).
	  gpInfo - the active set and other information stored in a
	   structure.
	  X - the input training data for the GP.
	  y - the output target training data for the GP.
		

	See also
	GPDECONSTRUCT, GPCREATE


	Copyright (c) 2007, 2009 Neil D. Lawrence
	
	"""
		
	options = gpOptions(gpInfo.approx);
	options.kern = kern;
	switch gpInfo.approx
	 case 'ftc'
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  options.numActive = size(gpInfo.X_u, 1);
	end
	model = gpCreate(size(X, 2), size(y, 2), X, y, options);
	model.scale = gpInfo.scale;
	model.bias = gpInfo.bias;
	model.m = gpComputeM(model);
	model.learnScales = gpInfo.learnScales;
	switch model.approx
	 case 'ftc'
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  model.beta = gpInfo.beta;
	  model.fixInducing = gpInfo.fixInducing;
	  if gpInfo.fixInducing
	    model.inducingIndices = gpInfo.inducingIndices;
	  else
	    model.X_u = gpInfo.X_u;
	  end
	end
	
	if gpInfo.d ~= size(y, 2)
	  error('y does not have correct number of dimensions.')
	end
	if gpInfo.q ~= size(X, 2)
	  error('X does not have correct number of dimensions.')
	end
	% FOrce update of everything.
	params = gpExtractParam(model);
	model = gpExpandParam(model, params);

def gpDisplay(model, spaceNum):

	"""Display a Gaussian process model.
	
	Description:
	
	gpDisplay(model, spaceNum) displays in human readable form the
	 contents of the GP model.
	 Arguments:
	  model - the model structure to be displaced.
	  spaceNum - number of spaces to place before displaying model
	   structure.
		

	See also
	GPCREATE, MODELDISPLAY.


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	if nargin > 1
	  spacing = repmat(32, 1, spaceNum);
	else
	  spaceNum = 0;
	  spacing = [];
	end
	spacing = char(spacing);
	fprintf(spacing);
	fprintf('Gaussian process model:\n')
	fprintf(spacing);
	fprintf('  Number of data points: %d\n', model.N);
	fprintf(spacing);
	fprintf('  Input dimension: %d\n', model.q);
	fprintf(spacing);
	fprintf('  Number of processes: %d\n', model.d);
	if isfield(model, 'beta') & ~isempty(model.beta)
	  fprintf(spacing);
	  fprintf('  beta: %2.4f\n', model.beta)
	end
	
	if any(model.scale~=1)
	  fprintf(spacing);
	  fprintf('  Output scales:\n');
	  for i = 1:length(model.scale)
	    fprintf(spacing);
	    fprintf('    Output scale %d: %2.4f\n', i, model.scale(i));
	  end
	end
	if any(model.bias~=0)
	  fprintf(spacing);
	  fprintf('  Output biases:\n');
	  for i = 1:length(model.bias)
	    fprintf(spacing);
	    fprintf('    Output bias %d: %2.4f\n', i, model.bias(i));
	  end
	end
	switch model.approx
	 case 'ftc'
	  fprintf(spacing);
	  fprintf('  No sparse approximation.\n')
	 case 'dtc'
	  fprintf(spacing);
	  fprintf('Deterministic training conditional approximation.\n')
	  fprintf('  Number of inducing variables: %d\n', model.k)
	 case 'dtcvar'
	  fprintf(spacing);
	  fprintf('Sparse variational approximation.\n')
	  fprintf('  Number of inducing variables: %d\n', model.k)
	 case 'fitc'
	  fprintf(spacing);
	  fprintf('Fully independent training conditional approximation.\n')
	  fprintf('  Number of inducing variables: %d\n', model.k)
	 case 'fitc'
	  fprintf(spacing);
	  fprintf('Partially independent training conditional approximation.\n')
	  fprintf('  Number of inducing variables: %d\n', model.k)
	end
	
	fprintf(spacing);
	fprintf('  Kernel:\n')
	kernDisplay(model.kern, 4+spaceNum);
	
	if isfield(model, 'noise') & ~isempty(model.noise)
	  fprintf(spacing);
	  fprintf('  Noise model:\n')
	  noiseDisplay(model.noise, 4+spaceNum);
	end
	

def modelRet = gpTest

	"""Test the gradients of the gpCovGrads function and the gp models.
	
	Description:
	
	model = gpTest runs some tests on the GP-LVM code in the GP
	 toolbox to test that it is working.
	 Returns:
	  model - a cell array of models used for testing.
		

	See also
	MODELTEST


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	q = 2;
	d = 3;
	N = 10;
	Nseq =4;
	k = 5;
	kernType = {'rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'};
	kernType = 'rbf';
	meanFunctionType = 'mlp';
	learnScales = true; % test learning of output scales.
	X = randn(N, q);
	Yorig = randn(N, d);
	indMissing = find(rand(N, d)>0.7);
	%indMissing = [9 19 29];
	approxType = {'ftc', 'dtc', 'dtcvar', 'fitc', 'pitc'};
	approxType = {'dtcvar'};
	counter = 0;
	for optimiseBeta = [false true]
	  for meanFunction = [false true]
	    %  for missing = [false true]
	    for missing = [false]
	      for fixInducing = [false true] 
	        Y = Yorig;
	        if missing
	          Y(indMissing) = NaN;
	        end
	        if meanFunction & missing
	          continue
	        end
	        for a = 1:length(approxType)
	          options = gpOptions(approxType{a});
	          options.learnScales = learnScales;
	          options.kern = kernType;
	          options.numActive = k;
	          options.isSpherical = ~missing;
	          options.isMissingData = missing;
	          options.fixInducing = fixInducing;
	          options.optimiseBeta = optimiseBeta;
	          if optimiseBeta & strcmp(approxType{a}, 'ftc')
	            options.beta = 1000;
	          end
	          if meanFunction 
	            disp(['Mean Function installed, ' ...
	                  'with ' approxType{a} ...
	                  ' approximation.'])
	          else
	            disp([approxType{a} ' approximation.'])  
	          end
	          if missing
	            disp(['Missing data used.'])
	          end
	          if fixInducing
	            disp(['Inducing variables fixed.'])
	            options.fixIndices = round(linspace(1, size(Y, 1), k));
	          end
	          if ~optimiseBeta
	            disp(['Beta not optimised.'])
	          end
	        
	          if meanFunction
	            options.meanFunction = meanFunctionType;
	            options.meanFunctionOptions = feval([meanFunctionType 'Options']);
	          end
	          model = gpCreate(q, d, X, Y, options);
	          
	          
	          initParams = gpExtractParam(model);
	          % this creates some nasty parameters.
	          initParams = randn(size(initParams));%./(100*randn(size(initParams)));
	          
	          % This forces kernel computation.
	          model = gpExpandParam(model, initParams);
	          gpCovGradsTest(model);
	          modelGradientCheck(model);
	          counter = counter + 1;
	          modelRet{counter} = model;
	        end
	      end
	    end
	  end
	end
def demGpCov2D(ind,bw):

	"""Simple demonstration of sampling from a covariance function.
	
	Description:
	
	demGpCov2D(index, bw) shows two dimensions of the covariance
	 matrix giving the joint distribution and the conditional.
	 DEMGPSAMPLE is run first to set the covariance matrix.
	 Arguments:
	  index - the indices of the two elements from the covariance matrix
	   for which the joint distribution should be displayed.
	  bw - whether or not to prepare the plot with a black and white
	   output (default false).

	See also
	DEMCOVFUNCSAMPLE, DEMGPSAMPLE


	Copyright (c) 2006, 2008 Neil D. Lawrence
	
	"""
		
	if nargin < 2
	  bw = false;
	end
	  
	global printDiagram
	if bw == true
	  contourColour = [0 0 0];
	  contourLineStyle = '--';
	  conditionalLineStyle = '-.';
	  conditionalLineColour = [0 0 0];
	  conditionalSize = 4;
	  fixedLineStyle = '-';
	  fixedLineColour = [0 0 0];
	  fixedSize = 4;
	else
	  contourColour = [0 0 1];
	  contourLineStyle = '-';
	  conditionalLineStyle = '-';
	  conditionalLineColour = [1 0 0];
	  conditionalSize = 4;
	  fixedLineStyle = '-';
	  fixedLineColour = [0 1 0];
	  fixedSize = 4;
	end
	contourSize = 4;
	if nargin < 1
	  ind = [1 2];
	end
	load demGpSample
	K = K(ind, ind);
	f = f(ind);
	x = x(ind);
	disp(K);
	figure(1)
	clf
	[ax, cont, t] = basePlot(K, ind);
	set(cont, 'linewidth', contourSize);
	set(cont, 'linestyle', contourLineStyle, 'color', contourColour);
	figure(2)
	clf
	[ax, cont, t] = basePlot(K, ind);
	set(cont, 'linewidth', contourSize);
	set(cont, 'linestyle', contourLineStyle, 'color', contourColour);
	
	cont2 = line([f(1) f(1)], [-1 1]);
	set(cont2, 'linewidth', fixedSize)
	set(cont2, 'linestyle', fixedLineStyle, 'color', fixedLineColour)
	
	figure(3)
	clf
	[ax, cont, t] = basePlot(K, ind);
	set(cont, 'linewidth', contourSize);
	set(cont, 'linestyle', contourLineStyle, 'color', contourColour);
	
	cont2 = line([f(1) f(1)], [-1 1]);
	set(cont2, 'linewidth', fixedSize)
	set(cont2, 'linestyle', fixedLineStyle, 'color', fixedLineColour)
	
	% Compute conditional mean and variance
	f2Mean = K(1, 2)/K(1, 1)*f(1);
	f2Var = K(2, 2) - K(1, 2)/K(1, 1)*K(1, 2);
	yval = linspace(-1, 1, 200);
	pdfVal = 1/sqrt(2*pi*f2Var)*exp(-0.5*(yval-f2Mean).*(yval-f2Mean)/f2Var);
	pdf = line(pdfVal*0.25, yval);
	set(pdf, 'linewidth', conditionalSize)
	set(pdf, 'linestyle', conditionalLineStyle, 'color', conditionalLineColour)
	
	if bw 
	  app = 'bw';
	else 
	  app = '';
	end
	for figNo = 1:3
	  figure(figNo)
	  if exist('printDiagram', 'var') & printDiagram
	    printPlot(['demGpCov2D' num2str(ind(1)) '_' num2str(ind(2)) '_' ...
	               num2str(figNo) app], '../tex/diagrams', '../html');
	  end
	end
	
def [ax, cont, t] = basePlot(K, ind):
	
	% BASEPLOT Plot the contour of the covariance.
	% FORMAT
	% DESC creates the basic plot.
	% 
	
	[U, V] = eig(K);
	r = sqrt(diag(V));
	theta = linspace(0, 2*pi, 200)';
	xy = [r(1)*sin(theta) r(2)*cos(theta)]*U;
	cont = line(xy(:, 1), xy(:, 2));
	
	
	t = text(0.5, -0.2, ['f_{' num2str(ind(1)) '}']);
	t =[t text(-0.2, -0.5, ['f_{' num2str(ind(2)) '}'])];
	
	
	set(gca, 'xtick', [-1  0  1])
	set(gca, 'ytick', [-1 0 1])
	
	set(t, 'fontname', 'times')
	set(t, 'fontsize', 24)
	set(t, 'fontangle', 'italic')
	
	zeroAxes(gca, 0.025, 18, 'times')
	
	ax = gca;
def [gmu, gsigmavar, factors] = gpPosteriorGradMeanCovar(model, X):

	"""Gadient of the mean and variances of the posterior at points given by X.
	
	Description:
	
	[gmu, gCovar] = gpPosteriorGradMeanCovar(model, X) computes the
	 gradient of the mean and covariances of the posterior distribution
	 of a Gaussian process with respect to the input locations.
	 Returns:
	  gmu - the gradient of the posterior mean with respect to the input
	   locations.
	  gCovar - the gradients of the posterior covariance with respect to
	   the input locations. By raw, we mean that the gradient has not yet
	   been multiplied by any output scale in that direction (as is done
	   for gpPosteriorGradMeanCovar). The gradients are stored in a cell
	   array of dimension MODEL.q x MODEL.d.
	 Arguments:
	  model - the model for which gradients are to be computed.
	  X - the input locations where gradients are to be computed.
		DESC computes the gradient of the mean and covariances of the
		posterior distribution of a Gaussian process with respect to the
		input locations. Returns a compact representation for the
		covariances which separates the factors associated with the
		different dimensions from the covariance gradients.
		ARG model : the model for which gradients are to be computed.
		ARG X : the input locations where gradients are to be computed.
		RETURN gmu : the gradient of the posterior mean with respect to
		the input locations.
		RETURN grCovar : the 'raw' gradient of the posterior covariance with
		respect to the input locations. By raw, we mean that the gradient
		has not yet been multiplied by any output scale in that direction
		(as is done for gpPosteriorGradMeanCovar). The gradients are
		stored in a cell array of length MODEL.q.
		RETURN factors : the factors for multiplying the 'raw' gradients
		of the covariances by.
		
		

	See also
	GPCREATE, GPPOSTERIORMEANVAR


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	if ~isfield(model, 'alpha')
	  model = gpComputeAlpha(model);
	end
	
	switch model.approx
	 case 'ftc'
	  gX = kernGradX(model.kern, X, model.X);
	  kX_star = kernCompute(model.kern, X, model.X)';
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  gX = kernGradX(model.kern, X, model.X_u);
	  kX_star = kernCompute(model.kern, X, model.X_u)';
	 otherwise
	  error('Unrecognised approximation type');
	  
	end
	K = kernGradX(model.kern, X);
	
	
	if ~model.isMissingData
	  for i = 1:model.q
	    switch model.approx
	     case 'ftc'
	      KinvgK = model.invK_uu*squeeze(gX(:, i, :));
	     case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	      KinvgK = (model.invK_uu - (1/model.beta)*model.Ainv)*squeeze(gX(:, i, :));
	     otherwise
	      error('Unrecognised approximation type');
	    end
	    kXTKinvgK = kX_star'*KinvgK;
	    gCovar{i} = squeeze(K(:, i, :))-kXTKinvgK - diag(diag(kXTKinvgK));
	    gmu{i} = squeeze(gX(:, i, :))'*model.alpha.*repmat(model.scale, ...
	                                                      size(X, 1), 1);
	  end
	  
	  
	  % Deal with scaling.
	  if nargout < 3
	    for i = 1:model.q
	      for j = 1:model.d
	        gsigmavar{i, j} = gCovar{i}*model.scale(j)*model.scale(j);
	      end
	    end
	  else
	    factors = model.scale.*model.scale;
	    gsigmavar = gCovar;
	  end
	else
	  error('Not yet implemented for models trained on missing data.');
	end

	"""Model silhouette data with independent linear models.
	
	Description:
	
	"""
		% FORMAT
	% DESC runs a simple regression on the Agawal and Triggs data.
	%
	% SEEALSO : demSilhouetteGp1, demSilhouetteAverage
	% 
	% COPYRIGHT : Neil D. Lawrence, 2008
	
	
	randn('seed', 1e7)
	rand('seed', 1e7)
	
	dataSetName = 'silhouette';
	experimentNo = 1;
	
	% load data
	[X, y, XTest, yTest] = mapLoadData(dataSetName);
	
	
	% Set up the model
	options = linearOptions;
	
	q = size(X, 2);
	d = size(y, 2);
	model = linearCreate(q, d, options);
	
	
	model = linearOptimise(model, X, y);
	modelDisplay(model)
	
	% Save results
	capName = dataSetName;;
	capName(1) = upper(capName(1));
	fileBaseName = ['dem' capName 'Linear' num2str(experimentNo)];
	save([fileBaseName '.mat'], 'model');
	
	
	demSilhouettePlot
def [mu, covarSigma, factors] = gpPosteriorMeanCovar(model, X):

	"""Mean and covariances of the posterior at points given by X.
	
	Description:
	
	[mu, Sigma] = gpPosteriorMeanCovar(model, X) gives the posterior
	 mean and covariance at the points given by X.
	 Returns:
	  mu - the posterior mean.
	  Sigma - the posterior covariance.
	 Arguments:
	  model - the model for which the posterior will be computed.
	  X - the latent positions where the mean and covariance will be
	   computed.
	
	[mu, Sigma, factor] = gpPosteriorMeanCovar(model, X) gives the
	 posterior mean and covariance at the points given by X without
	 scaling on the output posterior covariance. This allows for a more
	 compact representation. The scale factors are provided in a
	 separate vector FACTORS.
	 Returns:
	  mu - the posterior mean.
	  Sigma - the posterior covariance *without scaling*.
	  factor - the factors to multiply each dimension by to obtain the
	   covariances for each output.
	 Arguments:
	  model - the model for which the posterior will be computed.
	  X - the latent positions where the mean and covariance will be
	   computed.
		

	See also
	GPCREATE, GPPOSTERIORMEANVAR


	Copyright (c) 2006, 2009 Neil D. Lawrence
	
	"""
		
	% Just call gpPosteriorMeanVar for means.
	mu = gpPosteriorMeanVar(model, X);
	
	if nargout > 1
	  if size(X, 1)>1000
	    warning(['Computation of covariances takes a long time for larger ' ...
	             'data sets, are you sure you did''nt just want ' ...
	             'variances? If so use gpPosteriorMeanVar.'])
	  end
	end
	
	% Compute kernel for new point.
	switch model.approx
	 case 'ftc'
	  KX_star = kernCompute(model.kern, model.X, X);  
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  KX_star = kernCompute(model.kern, model.X_u, X);  
	end
	
	% Compute covariances if requried.
	if nargout > 1
	  % Compute kernel for new point.
	  K = kernCompute(model.kern, X);
	  if ~model.isMissingData
	    switch model.approx
	     case 'ftc'
	      Kinvk = model.invK_uu*KX_star;
	     case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	      Kinvk = (model.invK_uu - (1/model.beta)*model.Ainv)*KX_star;
	    end
	    
	    covarsig = K - KX_star'*Kinvk;
	    if isfield(model, 'beta')
	      covarsig = covarsig + eye(size(X, 1))*(1/model.beta);
	    end
	    if nargout>2
	      covarSigma = covarsig;
	      factors = model.scale.*model.scale;
	    else
	      for i = 1:model.d
	        covarSigma{i} = covarsig*model.scale(i)*model.scale(i);
	      end
	    end
	  else
	    for i = 1:model.d
	      switch model.approx
	       case 'ftc'
	        Kinvk{i} = model.invK_uu{i}*KX_star;
	       case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	        Kinvk{i} = (model.invK_uu - (1/model.beta)*model.Ainv{i})*KX_star;
	      end      
	      covarsig{i} = K - KX_star'*Kinvk{i};
	      if isfield(model, 'beta')
	        covarsig{i} = covarsig{i} + eye(size(X, 1))*(1/model.beta);
	      end
	      if nargout>2
	        error(['Cannot return covariance and scales with models ' ...
	               'trained on missing data.'])
	      else
	        covarSigma{i} = covarsig{i}*model.scale(i)*model.scale(i);
	      end
	    end
	  end
	end 

def ind = gpBlockIndices(model, blockNo):

	"""Return indices of given block.
	
	Description:
	
	indices = gpBlockIndices(model, blockNum) returns the indices of a
	 given block for the PITC approximation.
	 Returns:
	  indices - the data indices associated with given block.
	 Arguments:
	  model - the model for which the indices are being computed.
	  blockNum - the block number for which the indices are required.
		gpLogLikelihood, gpUpdateAD
		

	See also
	GPCOMPUTEALPHA, GPCOVGRADS, GPLOGLIKEGRADIENTS, 


	Copyright (c) 2006 Neil D. Lawrence
	
	"""
		
	
	if ~strcmp(model.approx, 'pitc')
	  error('Block indices only relevant for pitc approximation.');
	else
	  if blockNo == 1
	    startVal = 1;
	  else
	    startVal = model.blockEnd(blockNo-1)+1;
	  end
	  endVal = model.blockEnd(blockNo);
	  ind = startVal:endVal;
	end

def model = gpReadFromFID(FID, varargin):

	"""Load from a FID produced by the C++ implementation.
	
	Description:
	
	model = gpReadFromFID(FID) loads in from a file stream the data
	 format produced by the C++ GP implementation.
	 Returns:
	  model - the model loaded in from the file.
	 Arguments:
	  FID - the file ID from where the data is loaded.
		

	See also
	GPREADFROMFILE


	Copyright (c) 2005, 2006, 2008, 2009 Neil D. Lawrence
	
	"""
		
	numData = readIntFromFID(FID, 'numData');
	dataDim = readIntFromFID(FID, 'outputDim');
	inputDim = readIntFromFID(FID, 'inputDim');
	sparseApprox = readIntFromFID(FID, 'sparseApproximation');
	numActive = readIntFromFID(FID, 'numActive');
	if sparseApprox
	  beta = modelReadFromFID(FID);
	  beta = beta(1, 1);
	end
	
	learnScale = readBoolFromFID(FID, 'learnScale');
	learnBias = readBoolFromFID(FID, 'learnBias');
	scale = modelReadFromFID(FID);
	bias = modelReadFromFID(FID);
	
	kern = modelReadFromFID(FID);
	noise = modelReadFromFID(FID);
	
	if sparseApprox
	  X_u = modelReadFromFID(FID);
	end
	
	X = varargin{1};
	y = varargin{2};
	% X = zeros(numData, inputDim);
	% y = zeros(numData, dataDim);
	
	warning('Noise model is ignored');
	
	switch sparseApprox
	 case 0
	  approxType = 'ftc';
	 case 1
	  approxType = 'dtc';
	 case 2
	  approxType = 'fitc';
	 case 3
	  approxType = 'pitc';
	 case 4
	  approxType = 'dtcvar';
	end
	options = gpOptions(approxType);
	options.numActive = numActive;
	options.kern = kern;
	model = gpCreate(inputDim, size(y, 2), X, y, options);
	model.X = X;
	model.X_u = X_u;
	if sparseApprox
	  model.beta = beta;
	end
	model.scale = scale;
	model.bias = bias;
	model.m = gpComputeM(model);
	
	% This forces kernel computation.
	initParams = gpExtractParam(model);
	model = gpExpandParam(model, initParams);

def [mu, varsigma] = gpPosteriorMeanVar(model, X):

	"""Mean and variances of the posterior at points given by X.
	
	Description:
	
	[mu, sigma] = gpPosteriorMeanVar(model, x) returns the posterior
	 mean and variance for a given set of points.
	 Returns:
	  mu - the mean of the posterior distribution.
	  sigma - the variances of the posterior distributions.
	 Arguments:
	  model - the model for which the posterior will be computed.
	  x - the input positions for which the posterior will be computed.
		

	See also
	GPCREATE, GPPOSTERIORMEANCOVAR, GPPOSTERIORGRADMEANVAR


	Copyright (c) 2005, 2006, 2009 Neil D. Lawrence
	
	"""
		
	
	if ~isfield(model, 'alpha')
	  model = gpComputeAlpha(model);
	end
	
	maxMemory = 1000000;
	switch model.approx
	 case 'ftc'
	  chunkSize = ceil(maxMemory/model.N);
	 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	  chunkSize = ceil(maxMemory/model.k);
	end
	
	mu = zeros(size(X, 1), model.d);
	if nargout > 1
	  varsigma = zeros(size(X, 1), model.d);
	end
	
	startVal = 1;
	endVal = chunkSize;
	if endVal>size(X, 1)
	  endVal = size(X, 1);
	end
	
	while startVal <= size(X, 1)
	  indices = startVal:endVal;
	
	  % Compute kernel for new point.
	  switch model.approx
	   case 'ftc'
	    KX_star = kernCompute(model.kern, model.X, X(indices, :));  
	   case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	    KX_star = kernCompute(model.kern, model.X_u, X(indices, :));  
	  end
	
	  % Compute mean, using precomputed alpha vector.
	  if ~isfield(model, 'isMissingData') || ~model.isMissingData ... 
	    | ~strcmp(model.approx, 'ftc')
	    mu(indices, :) = KX_star'*model.alpha;
	  else
	    for i = 1:model.d
	      mu(indices, i) = KX_star(model.indexPresent{i}, :)' ...
	          *model.alpha(model.indexPresent{i}, i);
	    end
	  end
	  
	  % Compute variances if requried.
	  if nargout > 1
	    if ~isfield(model, 'isSpherical') || model.isSpherical
	      % Compute diagonal of kernel for new point.
	      diagK = kernDiagCompute(model.kern, X(indices, :));
	      switch model.approx
	       case 'ftc'
	        Kinvk = model.invK_uu*KX_star;
	       case {'dtc', 'dtcvar', 'fitc', 'pitc'}
	        Kinvk = (model.invK_uu - (1/model.beta)*model.Ainv)*KX_star;
	      end
	      varsig = diagK - sum(KX_star.*Kinvk, 1)';
	      if isfield(model, 'beta')
	        varsig = varsig + (1/model.beta);
	      end
	      varsigma(indices, :) = repmat(varsig, 1, model.d);
	    else
	      diagK = kernDiagCompute(model.kern, X(indices, :));
	      for i = 1:model.d
	        ind = model.indexPresent{i};
	        switch model.approx
	         case 'ftc'
	          Kinvk = model.invK_uu{i}*KX_star(ind, :);
	         otherwise 
	          error(['Non-spherical not yet implemented for any approximation ' ...
	                 'other than ''ftc''']);
	        end
	        varsigma(indices, i) = diagK - sum(KX_star(ind, :).*Kinvk, 1)';
	      end
	    end
	  end
	  
	  
	  % Rescale the mean
	  mu(indices,:) = mu(indices, :).*repmat(model.scale, length(indices), 1);
	  % Add the bias back in.
	  mu(indices, :) = mu(indices, :) + repmat(model.bias, length(indices), 1);
	  % If the mean function is present, add it it.
	  if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
	    mu(indices, :) = mu(indices, :) + modelOut(model.meanFunction, ...
	                                               X(indices, :));
	  end
	  % rescale the variances
	  if nargout > 1
	    varsigma(indices, :) = varsigma(indices, :).* ...
	        repmat(model.scale.*model.scale, length(indices), 1);
	  end
	
	  % Prepare for the next chunk.
	  startVal = endVal + 1;
	  endVal = endVal + chunkSize;
	  if endVal > size(X, 1)
	    endVal = size(X, 1);
	  end
	end
def g = gpMeanFunctionGradient(model):

	"""Compute the log likelihood gradient wrt the scales.
	
	Description:
	
	gpMeanFunctionGradient(model, g) computes the gradient of the log
	 likelihood with respect to the scales. In the future the gradients
	 with respect to the biases may also be included.
	 Arguments:
	  model - the model for which the gradients are to be computed.
	  g - the gradients of the likelihood with respect to the mean
	   function's parameters.
		

	See also
	GPCREATE, GPSCALEBIASGRADIENT, GPLOGLIKEGRADIENTS, GPLOGLIKELIHOOD


	Copyright (c) 2006, 2009 Neil D. Lawrence
	
	"""
		
	if isfield(model, 'isSpherical') & ~model.isSpherical
	  error('Currently only implemented for spherical');
	else
	  if model.isMissingData
	    error('Currently not implemented for missing data.');
	  end
	end
	  
	if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
	  g = zeros(1, model.meanFunction.numParams);
	  % compute gradients here.
	  switch model.approx
	   case 'ftc'
	    gmu = model.invK_uu*model.m;
	   case {'dtc', 'dtcvar'}
	    gmu = (model.m - model.K_uf'*model.Ainv*(model.K_uf*model.m))*model.beta;
	   case 'fitc'
	    Dinvm = model.Dinv*model.m;
	    gmu = (Dinvm-(model.Dinv*model.K_uf')...
	           *(model.Ainv*model.K_uf)*Dinvm)*model.beta;   
	   case 'pitc'
	    % Loop through the blocks computing each part to be added.
	    gmu = zeros(model.N, model.d);
	    K_ufDinvm = zeros(model.k, model.d);
	    K_ufDinv = zeros(model.k, model.N);
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      Dinvm{i} = model.Dinv{i}*model.m(ind, :);
	      K_ufDinvm = K_ufDinvm + model.K_uf(:, ind)*Dinvm{i};
	    end
	    for i = 1:length(model.blockEnd)
	      ind = gpBlockIndices(model, i);
	      gmu(ind, :) = (Dinvm{i} ...
	                     - model.Dinv{i}...
	                     *model.K_uf(:, ind)'*(model.Ainv*K_ufDinvm))*model.beta;
	    end
	  end
	  gmu = gmu./repmat(model.scale, model.N, 1);
	  goutputDparam = modelOutputGrad(model.meanFunction, model.X);
	  for i = 1:model.meanFunction.numParams
	    g(1, i) = sum(sum(gmu.*squeeze(goutputDparam(:, i, :))));
	  end
	else
	  g = [];
	end

def model = gpOptimise(model, display, iters,gradcheck):

	"""Optimise the inducing variable based kernel.
	
	Description:
	
	model = gpOptimise(model, display, iters) optimises the Gaussian
	 process  model for a given number of iterations.
	 Returns:
	  model - the optimised model.
	 Arguments:
	  model - the model to be optimised.
	  display - whether or not to display while optimisation proceeds,
	   set to 2 for the most verbose and 0 for the least verbose.
	  iters - number of iterations for the optimisation.
		
		

	See also
	SCG, CONJGRAD, GPCREATE, GPGRADIENT, GPOBJECTIVE


	Copyright (c) 2005, 2006 Neil D. Lawrence
	
	
	With modifications by Carl Henrik Ek 2008
	
	"""
		
	if(nargin<4)
	  gradcheck = false;
	  if nargin < 3
	    iters = 2000;
	    if nargin < 2
	      display = 1;
	    end
	  end
	end
	
	
	params = gpExtractParam(model);
	
	options = optOptions;
	if display
	  options(1) = 1;
	  if length(params) <= 100 && gradcheck
	    options(9) = 1;
	  end
	end
	options(14) = iters;
	
	if isfield(model, 'optimiser')
	  optim = str2func(model.optimiser);
	else
	  optim = str2func('conjgrad');
	end
	
	if strcmp(func2str(optim), 'optimiMinimize')
	  % Carl Rasmussen's minimize function 
	  params = optim('gpObjectiveGradient', params, options, model);
	else
	  % NETLAB style optimization.
	  params = optim('gpObjective', params,  options, ...
	                 'gpGradient', model);
	end
	
	model = gpExpandParam(model, params);

