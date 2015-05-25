def X = doubleMatrixReadFromFID(FID):

    """Read a full matrix from an FID.
    
    Description:
    
    X = doubleMatrixReadFromFID(FID) reads a matrix from an FID.
     Returns:
      X - the returned matrix read from the file.
     Arguments:
      FID - the file ID to read the matrix from.
        

    See also
    MODELREADFROMFILE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    numRows = readIntFromFID(FID, 'numRows');
    numCols = readIntFromFID(FID, 'numCols');
    
    X = zeros(numRows, numCols);
    for i = 1:numRows
      lineStr = getline(FID);
      tokens = tokenise(lineStr, ' ');
      if strcmp(tokens(end),'');
        tokens=tokens(1:end-1);
      end
      if(length(tokens)~=numCols)
        error('Incorrect file format.');
      end
      for j = 1:numCols
        X(i,j) = str2num(tokens{j});
      end
    end
def doubleMatrixWriteToFID(val, FID):

    """Writes a double matrix to an FID.
    
    Description:
    
    doubleMatrixWriteToFID(val, FID) writes a double matrix to a
     stream.
     Arguments:
      val - matrix to place in file.
      FID - stream to write to.
        

    See also
    DOUBLEMATRIXREADFROMFID, MATRIXWRITETOFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    writeVersionToFID(FID, 0.2);
    writeStringToFID(FID, 'baseType', 'matrix');
    writeStringToFID(FID, 'type', 'doubleMatrix');
    writeIntToFID(FID, 'numRows', size(val, 1));
    writeIntToFID(FID, 'numCols', size(val, 2));
    for i = 1:size(val, 1)
      for j = 1:size(val, 2)-1
        fprintf(FID, '%1.17e ', val(i, j));
      end
      fprintf(FID, '%1.17e\n', val(i, end));
    end
      

def neighboursInd = findNeighbours(Y, k):

    """find the k nearest neighbours for each point in Y.
    
    Description:
    
    findNeighbours(y, k) returns the indices of the k nearest
     neighbours to each point in the given data matrix Y.
     Arguments:
      y - the data in which neighbours need to be found.
      k - the number of neighbours that need to be found.
        

    See also
    LLECREATE, MVUCREATE, ISOMAPCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    Y2 = sum(Y.*Y, 2);
    D = repmat(Y2', size(Y, 1), 1) + repmat(Y2, 1, size(Y, 1)) - 2*Y*Y';
    D(1:size(D, 1)+1:end) = inf;
    [void, ind] = sort(D);
    neighboursInd = ind(1:k, :)';

def handle = imageModify(handle, imageValues, imageSize, transpose, negative, ...
		     scale):

    """Helper code for visualisation of image data.
    
    Description:
    
    handle = imageModify(handle, imageValues, imageSize, transpose,
     negative, scale) is a helper function for visualising image data
     using latent variable models.
     Returns:
      handle - a the handle to the image data.
     Arguments:
      handle - the handle of the image data.
      imageValues - the values to set the image data to.
      imageSize - the size of the image.
      transpose - whether the resized image needs to be transposed
       (default 1, which is yes).
      negative - whether to display the negative of the image (default
       0, which is no).
      scale - dummy input, to maintain compatability with
       IMAGEVISUALISE.
        

    See also
    IMAGEVISUALISE, FGPLVMRESULTSDYNAMIC


    Copyright (c) 2003, 2004, 2006 Neil D. Lawrence
    
    """
        
    
    if nargin < 4
      transpose = 1;
    end
    if nargin< 5
      negative = 0;
    end
    if negative
      imageValues = -imageValues;
    end
    if transpose
      set(handle, 'CData', reshape(imageValues, imageSize(1), imageSize(2))');
    else
      set(handle, 'CData', reshape(imageValues, imageSize(1), imageSize(2)));
    end

def handle = imageVisualise(imageVals, imageSize, transpose, negative, ...
				 scale):

    """Helper code for showing an image during 2-D visualisation.
    
    Description:
    
    handle = imageVisualise(imageValues, imageSize, transpose,
     negative, scale) is a helper function for plotting image data
     using latent variable models.
     Returns:
      handle - a the handle to the image data.
     Arguments:
      imageValues - the values to set the image data to.
      imageSize - the size of the image.
      transpose - whether the resized image needs to be transposed
       (default 1, which is yes).
      negative - whether to display the negative of the image (default
       0, which is no).
      scale - whether or not to use the imagesc function (defaults to 1,
       which is yes).
        

    See also
    IMAGEMODIFY, FGPLVMRESULTSDYNAMIC


    Copyright (c) 2003, 2004, 2006 Neil D. Lawrence
    
    """
        
    if nargin < 3
      transpose = 1;
    end
    if nargin< 4
      negative = 0;
    end
    if nargin < 5
      scale = 1;
    end
    if negative
      imageVals = -imageVals;
    end
    imageData = reshape(imageVals, imageSize(1), imageSize(2));
    if transpose
      imageData = imageData';
    end
    if scale
      handle = imagesc(imageData);
    else
      handle = image(imageData);
    end
def [X, sigma2] = isomapEmbed(Y, dims):

    """Embed data set with Isomap.
    
    Description:
    def [X, sigma2] = isomapEmbed(Y, dims):
%
    """
        
    % Note: isomap code uses the transpose of a design matrix.
    if any(any(isnan(Y)))
      error('Cannot initialise gplvm using isomap when missing data is present.')
    else
      D = L2_distance(Y', Y', 1);
      options.dims = 1:dims;
      neighbours = 7;
      [Xstruct, sigma2, E] = Isomap(D, 'k', neighbours, options);
      X = zeros(size(Y, 1), 2);
      if length(Xstruct.index) ~= size(Y, 1)
        % We don't really deal with this problem correctly here ...
        warning('Isomap graph is not fully connected');
      end
      X(Xstruct.index, :) = Xstruct.coords{dims}';
      % Rescale X so that variance is 1 and mean is zero.
      meanX = mean(X);
      X = X-ones(size(Y, 1), 1)*meanX;
      varX = var(X);
      X = X*diag(sqrt(1./varX));
    end
def model = kbrCreate(inputDim, outputDim, options):

    """Create a KBR model.
    
    Description:
        The kernel based regression model is simply a model for least
        squares regression in a kernel feature space. Any kernel from the KERN
        toolbox can be specified. The model was developed for providing kernel
        based back constraints in the GP-LVM. Please consider using a Gaussian
        process model (through the GP toolbox) if you are interested in the
        model for regression.
        
    
    model = kbrCreate(options) creates a kernel based regression model
     structure given an options structure.
     Returns:
      model - the model structure with the default parameters placed in.
     Arguments:
      options - an options structure that determines the form of the
       model.
        

    See also
    KBROPTIONS, KBRPARAMINIT, MODELCREATE


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    
    model.type = 'kbr';
    model.inputDim = inputDim;
    model.outputDim = outputDim;
    model.numData = size(options.X, 1);
    model.numParams = (model.numData + 1)*outputDim;
    model.X = options.X;
    if isstruct(options.kern)
      model.kern = options.kern;
    else
      model.kern = kernCreate(options.X, options.kern);
    end
    
    model.K = kernCompute(model.kern, options.X);
    model = kbrParamInit(model);

def kbrDisplay(model, spacing):

    """Display parameters of the KBR model.
    
    Description:
    
    kbrDisplay(model) displays the parameters of the kernel based
     regression model and the model type to the console.
     Arguments:
      model - the model to display.
    
    kbrDisplay(model, spacing)
     Arguments:
      model - the model to display.
      spacing - how many spaces to indent the display of the model by.
        

    See also
    KBRCREATE, MODELDISPLAY


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Kernel based regression model:\n')
    fprintf(spacing);
    fprintf('Kernel type:\n')
    kernDisplay(model.kern, length(spacing)+2)

def model = kbrExpandParam(model, params,dim):

    """Create model structure from KBR model's parameters.
    
    Description:
    
    model = kbrExpandParam(model, params) takes a vector of KBR
     weights and centres and places them in their respective positions
     in the KBR model.
     Returns:
      model - the model with the weights distributed in the correct
       places.
     Arguments:
      model - the model in which the weights are to be placed.
      params - a vector of the weights to be placed in the model.
        
        
        

    See also
    KBRUNPAK, KBRCREATE, KBREXTRACTPARAM


    Copyright (c) 2008 Neil D. Lawrence
    
    
    With modifications by Carl Henrik Ek 2008
    
    """
        if(nargin<3)
      startVal = 1;
      endVal = model.numData*model.outputDim;
      model.A = reshape(params(1:endVal), model.numData, model.outputDim);
      model.bias = params(endVal+1:end);
    else
      model.A(:,dim) = params(1:1:model.numData*length(dim));
      model.bias(dim) = params(model.numData*length(dim)+1:1:end);
    end
def [params, names] = kbrExtractParam(model,dim):

    """Extract parameters from the KBR model structure.
    
    Description:
    
    param = kbrExtractParam(model) extracts parameters from the kernel
     based regression model structure into a vector of parameters for
     optimisation.
     Returns:
      param - vector of parameters extracted from the model.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
        DESC extracts parameters and parameter names from the kernel based regression model structure.
        ARG model : the model structure containing the parameters to be
        extracted.
        RETURN param : vector of parameters extracted from the model.
        RETURN names : cell array of strings containing names for each parameter.
        
        
        
        
        

    See also
    KBRCREATE, KBREXPANDPARAM, MODELEXTRACTPARAM, SCG, CONJGRAD


    Copyright (c) 2005, 2006, 2008 Neil D. Lawrence
    
    
    With modifications by Carl Henrik Ek 2007
    
    """
        
    if(nargin<2)
      params = [model.A(:)' model.bias];
    else
      params = model.A(:,dim);
      params = [params(:)' model.bias(dim)];
    end
    
    if nargout > 1
      % Add names to parameters
      counter = 0;
      for i = 1:model.numData
        for j = 1:model.outputDim
          counter = counter + 1;
          names{counter} = ['A(' num2str(i) ', ' num2str(j) ')'];
        end
      end
      for j = 1:model.outputDim
        counter = counter + 1;
        names{counter} = ['bias(' num2str(j) ')'];
      end
    end
def model = kbrOptimise(model, X, Y, varargin):

    """Optimise a KBR model.
    
    Description:
    
    model = kbrOptimise(model, X, Y) optimises a kernel based
     regression model using a least squares fit.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
      X - the input data locations for the optimisation.
      Y - the target data locations for the optimisation.
        

    See also
    KBRCREATE, MODELOPTIMISE


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    
    model.numData = size(X, 1);
    model.K = kernCompute(model.kern, X);
    model.X = X;
    
    model.bias = mean(Y, 1);
    model.A = pdinv(model.K)*(Y-repmat(model.bias, model.numData, 1));

def options = kbrOptions(X):

    """Create a default options structure for the KBR model.
    
    Description:
    
    options = kbrOptions(X) creates a default options structure for
     the kernel based regression model structure.
     Returns:
      options - the default options structure.
     Arguments:
      X - the input data for the kernel regression.
        

    See also
    KBRCREATE, MODELOPTIONS


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    options.kern = 'rbf';
    options.X = X;

def [Y, G] = kbrOut(model, X):

    """Compute the output of a KBR model given the structure and input X.
    
    Description:
    
    y = kbrOut(model, x) computes the model parameters for the kernel
     based regression model given inputs associated with rows and
     columns.
     Returns:
      y - the output results.
     Arguments:
      model - the model structure for which the output is computed.
      x - the input data.
    
    [Y, G] = kbrOut(model, X) gives the output of a radial basis
     function model.
     Returns:
      Y - the output.
      G - the values computed at the kernel.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
        

    See also
    KBRCREATE, MODELCOMPUTE, MODELCREATE, KBREXPANDPARAM, KBREXTRACTPARAM


    Copyright (c) 2005, 2006, 2008 Neil D. Lawrence
    
    """
        
    numData = size(X, 1);
    if ~isfield(model, 'bias') & isfield(model, 'b')
      model.bias = model.b;
    end
    G = kernCompute(model.kern, X, model.X);
    Y = G*model.A+ones(numData, 1)*model.bias;

def g = kbrOutputGrad(model, X, dim):

    """Evaluate derivatives of KBR model outputs with respect to parameters.
    
    Description:
    
    g = kbrOutputGrad(model, X) evaluates the derivates of a kernel
     based regression outputs with respect to the parameters of the
     multi-layer perceptron.
     Returns:
      g - the gradient of the outputs of the kernel based regression
       with respect to each of the parameters. The size of the matrix is
       number of data x number of parameters x number of outputs of the
       model.
     Arguments:
      model - the model for which the derivatives are to be computed.
      X - the input data locations where the gradients are to be
       computed.
        
        

    See also
    KBRCREATE, KBRDERIV


    Copyright (c) 2008 Neil D. Lawrence
    
    
    With modifications by Carl Henrik Ek 2008
    
    """
        
    numData = size(X, 1);
    if(nargin<=2)
      for i = 1:model.outputDim
        startZeros = zeros(numData, numData*(i - 1));
        finishZeros = zeros(numData, numData*(model.outputDim-i));
        startZeros2 = zeros(numData, (i - 1));
        finishZeros2 = zeros(numData, (model.outputDim-i));
        g(:, :, i) = [startZeros model.K finishZeros startZeros2 ones(numData, 1) finishZeros2];
      end
    else
      g(:,:) = [model.K ones(numData,1)];
    end
def model = kbrParamInit(model):

    """KBR model parameter initialisation.
    
    Description:
    
    model = kbrParamInit(model) initialises the kernel based
     regression model structure with some default parameters.
     Returns:
      model - the model structure with the default parameters placed in.
     Arguments:
      model - the model structure which requires initialisation.
        

    See also
    KBRCREATE, MODELCREATE, MODELPARAMINIT


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    model.A = randn(model.numData, model.outputDim)/sqrt(model.numData+1);
    model.bias = randn(1, model.outputDim)/sqrt(model.numData+1);
    
    
    

def [X, sigma2] = kpcaEmbed(Y, dims):

    """Embed data set with kernel PCA.
    
    Description:
    def [X, sigma2] = kpcaEmbed(Y, dims):
%
    """
        
    
    if any(any(isnan(Y)))
      error('When missing data is present Kernel PCA cannot be used to initialise')
    end
    
    K = kernCompute(kern, Y);
    [u, v] = eigs(K, dims);
    X = u*sqrt(v);
    sigma2 = -1;
def lfmClassVisualise( call ):

    """Callback function to visualize LFM in 2D
    
    Description:
    def lfmClassVisualise( call ):
%
    """
        
    global visualiseInfo
    maxN = 30;
    
    switch call
     case 'click'       
      if ~visualiseInfo.clicked  
        visualiseInfo.clicked = ~visualiseInfo.clicked;       
        tic;
        visualiseInfo.lastToc = 0;
        visualiseInfo.timer.series = [];
        visualiseInfo.f1.series = [];
        visualiseInfo.f2.series = [];
      else
        visualiseInfo.clicked = ~visualiseInfo.clicked;       
        N = length(visualiseInfo.timer.series);
        if N >maxN
          ind = round(linspace(1, N, maxN));
        else
          ind = 1:N;
        end
         %disp(ind)
         xVector = visualiseInfo.timer.series(ind)';
         f{1} = visualiseInfo.f1.series(ind)';
         f{2} = visualiseInfo.f2.series(ind)';
         %xVector = linspace(0,3,30)';
         %f{1} = visualiseInfo.varargin{4}{1};
         %f{2} = visualiseInfo.varargin{4}{2};     
         Y = modelOut(visualiseInfo.model, xVector, f, visualiseInfo.varargin{1});
         if strcmp(visualiseInfo.model.approx, 'ftc')
             channels = repmat(visualiseInfo.varargin{1}, length(ind), 1);
             channelsLabels = [41:47 49:50];
             for k = 1: length(channelsLabels)
                 channels(:, channelsLabels(k)) = Y(:,k);
             end
         else
             channels = Y;         
         end
         for j = 1:size(channels, 1)
           if j>1
           pause(xVector(j) -xVector(j-1));
           end
           visualiseInfo.visualiseModify(visualiseInfo.visHandle, channels(j, :) , ... 
                                         visualiseInfo.varargin{2});
         end
         visualiseInfo.timer.series = 0;
         visualiseInfo.f1.series = 0;
         visualiseInfo.f2.series = 0;
         %set(visualiseInfo.f1.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f1.series);
         %set(visualiseInfo.f2.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f2.series);     
         
        end                
        case 'move'
            if visualiseInfo.clicked
              timeNow = toc;
              if ~isfield(visualiseInfo, 'lastToc') || timeNow > visualiseInfo.lastToc + 1/24
                visualiseInfo.lastToc = timeNow;
              [x, y]  = localCheckPointPosition(visualiseInfo);
                if ~isempty(x)              
                    set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);                               
                    visualiseInfo.timer.series = [visualiseInfo.timer.series timeNow];
                    visualiseInfo.f1.series = [visualiseInfo.f1.series x];           
                    visualiseInfo.f2.series = [visualiseInfo.f2.series y];           
                    %set(visualiseInfo.f1.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f1.series);               
                    %set(visualiseInfo.f2.handle,'Xdata', visualiseInfo.timer.series, 'Ydata', visualiseInfo.f2.series);               
                    %sprintf('%f\n',length(visualiseInfo.f1.series))                
                end            
            end
            end
        otherwise
    
        %
    end
    
def [x, y] = localCheckPointPosition(visualiseInfo):
    
    % Get the point of the cursor
    point = localGetNormCursorPoint(gcf);
    
    % get the position of the axes
    position = get(visualiseInfo.plotAxes, 'Position');
    
    % Check if the pointer is in the axes
    if point(1) > position(1) ...
          && point(1) < position(1) + position(3) ...
          && point(2) > position(2) ...
          && point(2) < position(2) + position(4);
      
      % Rescale the point according to the axes
      [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);
    
      % Find the nearest point
    else
      % Return nothing
      x = [];
      y = [];
    end
    
def point = localGetNormCursorPoint(figHandle):
    
    point = get(figHandle, 'currentPoint');
    figPos = get(figHandle, 'Position');
    % Normalise the point of the curstor
    point(1) = point(1)/figPos(3);
    point(2) = point(2)/figPos(4);
    
def [x, y] = localGetNormAxesPoint(point, axesHandle):
    
    position = get(axesHandle, 'Position');
    x = (point(1) - position(1))/position(3);
    y = (point(2) - position(2))/position(4);
    lim = get(axesHandle, 'XLim');
    x = x*(lim(2) - lim(1));
    x = x + lim(1);
    lim = get(axesHandle, 'YLim');
    y = y*(lim(2) - lim(1));
    y = y + lim(1);
    
    

def lfmResultsDynamic(dataSet, number, dataType, varargin):

    """Load a results file and visualise them.
    
    Description:
    
    lfmResultsDynamic(dataSet, number, dataType, ...) loads results of
     a latent variable model and visualises them.
     Arguments:
      dataSet - the name of the data set to load.
      number - the number of the run used.
      dataType - the type of data to visualise.
      ... - additional arguments to be passed to the lvmVisualise
       command.
        
        

    See also
    LFMVISUALISE


    Copyright (c) 2008 Neil D. Lawrence
    
    
    With modifications by Mauricio Alvarez 2009
    
    """
        
    %[model, lbls] = lvmLoadResult(modelType, dataSet, number);
    capName = dataSet;
    capName(1) = upper(capName(1));
    load(['dem' capName num2str(number) '.mat'], 'model');
    %load(['additional' capName '.mat'], 'skel', 'initPos');  
    
    % Visualise the results
    switch model.nlf
     case 2
      lfmVisualise(model, [dataType 'Visualise'], [dataType 'Modify'], varargin{:});
      
     otherwise 
      error('No visualisation code for data with this number of latent forces.');
    end
def lfmVisualise(model, visualiseFunction, visualiseModify, varargin):

    """Visualise the outputs in a latent force model
    
    Description:
    
    lfmVisualise(model, visualiseFunction, visualiseModify, ...)
     visualises a latent force model with two latent forces as inputs
     Arguments:
      model - the model to visualise.
      visualiseFunction - the function that draws the visualisation (in
       data space) when the graphs are first drawn.
      visualiseModify - the function that modifies the visualisation as
       you create the latent function.
      ... - various additional arguments to be passed to the
       visualisation commands.
        


    Copyright (c) 2008 Mauricio Alvarez and Neil D. Lawrence
    
    """
        
    global visualiseInfo
    kernType = model.kernType(1:3);
    if ~strcmpi(kernType, 'lfm')
        error('This function is only implemented for "LFM" kernels');
    else
        if length(varargin)~=4,
           error('Include the skeleton, the initial position and two posteriors for the outputs')
        end     
    end
    range = 1.1*max(max(abs([varargin{3}{1} varargin{3}{2}])));
    figure(1)
    set(gcf,'Position',[19 202 582 527]);
    clf
    % Create a black panel as background
    f1 = linspace(-range, range, 200);
    f2 = linspace(-range, range, 200);
    ax = axes('position', [0.15 0.15 0.75 0.75]);
    hold on
    sizeBack = 200;
    mapBack = repmat([0.9 0.9 0.9], sizeBack*sizeBack , 1);
    hold off
    f1Lim = [min(f1) max(f1)];
    f2Lim = [min(f2) max(f2)];
    set(ax, 'xLim', f1Lim);
    set(ax, 'yLim', f2Lim);
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 15);
    set(ax, 'TickLength', [0 0]);
    plot(varargin{3}{1}, varargin{3}{2}, 'xb', 'markersize', 10, 'LineWidth', 2)
    hold on
    plot(varargin{4}{1}, varargin{4}{2}, 'xr', 'markersize', 10, 'LineWidth', 2)
    xlabel('f_1(t)')
    ylabel('f_2(t)')
    grid on
    hold off
    %
    visualiseInfo.plotAxes = ax; 
    visualiseInfo.latentHandle = line(varargin{3}{1}(1), varargin{3}{2}(1), 'markersize', 20, 'color', ...
                                     [0 0 0], 'marker', '.', 'visible', ...
                                     'on', 'erasemode', 'xor');
    visualiseInfo.clicked = 0;
    visualiseInfo.digitAxes = [];
    visualiseInfo.digitIndex = [];
    
    visualiseInfo.dynamicsSlider = ...
        uicontrol('Style', 'slider', ...
                  'String', 'Time', ...
                  'sliderStep', [0.01, 0.1], ...
                  'units', 'normalized', ...
                  'position', [0 0.95 1 0.05], ...
                  'callback', 'lfmClassVisualise(''dynamicsSliderChange'')');
    
    set(visualiseInfo.dynamicsSlider, 'visible', 'off');
    set(gcf, 'WindowButtonMotionFcn', 'lfmClassVisualise(''move'')')
    set(gcf, 'WindowButtonDownFcn', 'lfmClassVisualise(''click'')')
    
    visualiseInfo.timer.series = 0;
    visualiseInfo.timer.stepTime = 0.03;
    
    %figure(2)
    %set(gcf,'Position',[418 332 445 395]);
    %clf
    %subplot(2,1,1);
    %visualiseInfo.f1.handle = plot(0,1.5, 'LineWidth', 2);
    %title('f_1(t)','FontSize', 15, 'FontName', 'arial')
    %set(gca, 'yLim', [-range range]);
    %set(gca, 'fontname', 'arial');
    %set(gca, 'fontsize', 15);
    visualiseInfo.f1.series = 0;
    %subplot(2,1,2);
    %visualiseInfo.f2.handle = plot(0,1.5,'LineWidth', 2);
    %title('f_2(t)','FontSize', 15, 'FontName', 'arial')
    %set(gca, 'yLim', [-range range]);
    %set(gca, 'fontname', 'arial');
    %set(gca, 'fontsize', 15);
    visualiseInfo.f2.series = 0;
    
    figure(2)
    set(gcf,'Position',[635 201 582 527]);
    clf
    visualiseInfo.visualiseFunction = str2func(visualiseFunction);
    visHandle = visualiseInfo.visualiseFunction(varargin{1:2});
    %set(gca, 'xlim', [-8  18], ...
    %         'ylim', [-2 15], ...
    %         'zlim', [0 35]);
    set(gca, 'xlim', [-20 10], ...
             'ylim', [-20 20], ...
             'zlim', [0 32]);
    title('Synthetic Motion', 'FontSize', 15);
         
    % Pass the data to visualiseInfo
    visualiseInfo.model = model;
    visualiseInfo.varargin = varargin;
    visualiseInfo.visualiseModify = str2func(visualiseModify);
    visualiseInfo.visHandle = visHandle;
    
    % figure(4)
    % clf
    % set(gcf,'Position',[418 30 445 395]);
    % visualiseInfo.fPos = plot(0,1.5,'LineWidth', 2);
    % %set(gca, 'xLim', [-2*range 2*range]);
    % %set(gca, 'yLim', [-2*range 2*range]);
    % xlabel('f_1(t)')
    % ylabel('f_2(t)')
    % set(gca, 'fontname', 'arial');
    % set(gca, 'fontsize', 15);
    
                                  
def model = linearCreate(inputDim, outputDim, options):

    """Create a linear model.
    
    Description:
    def model = linearCreate(inputDim, outputDim, options):
%
    """
        
    model.type = 'linear';
    model.activeFunc = options.activeFunc; 
    model.inputDim = inputDim;
    model.outputDim = outputDim;
    model.numParams = (inputDim + 1)*outputDim;
    
    model = linearParamInit(model);

def linearDisplay(model, spacing):

    """Display a linear model.
    
    Description:
    def linearDisplay(model, spacing):
%
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Model model:\n')
    fprintf(spacing);
    fprintf('  Input dimension: %d\n', model.inputDim);
    fprintf(spacing);
    fprintf('  Output dimension: %d\n', model.outputDim);

def model = linearExpandParam(model, params):

    """Update linear model with vector of parameters.
    
    Description:
    def model = linearExpandParam(model, params):
%
    """
        
    
    startVal = 1;
    endVal = model.inputDim*model.outputDim;
    model.W = reshape(params(1:endVal), model.inputDim, model.outputDim);
    model.b = params(endVal+1:end);
def [params, names] = linearExtractParam(model):

    """Extract weights from a linear model.
    
    Description:
    def [params, names] = linearExtractParam(model):
%
    """
        
      params = [model.W(:)' model.b];
    if nargout > 1
      counter = 0;
      for j = 1:size(model.W, 2)
        for i = 1:size(model.W, 1)
          counter = counter + 1;
          names{counter} = ['Weight ' num2str(i) '-' num2str(j)];
        end
      end
        for j = 1:size(model.b, 2)
        counter = counter + 1;
        names{counter} = ['Bias ' num2str(j)];
      end
    end
      
def g = linearLogLikeGradients(model):

    """Linear model gradients.
    
    Description:
    
    g = linearLogLikeGradients(model) computes the gradients of the
     log likelihood of a linear model with respect to the parameters.
     Returns:
      g - the gradients of the model log likelihood.
     Arguments:
      model - the model structure for computing the log likelihood.
        

    See also
    MODELLOGLIKEIHOOD, LINEARGRAD


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    Xo = [model.X ones(size(model.X, 1), 1)];
    W = [model.W; model.b];
    G = -model.beta*(Xo'*Xo*W - Xo'*model.y);
    gW = G(1:end-1, :);
    gb = G(end, :);
    g = [gW(:)' gb];

def ll = linearLogLikelihood(model):

    """Linear model log likelihood.
    
    Description:
    
    ll = linearLogLikelihood(model) computes the log likelihood of a
     linear model.
     Returns:
      ll - the model log likelihood.
     Arguments:
      model - the model structure for computing the log likelihood.
        

    See also
    MODELLOGLIKEIHOOD, MLPERR


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    N = size(model.y, 1);
    centred = model.y - repmat(model.b, N, 1) - model.X* ...
           model.W;
    
    centred = centred*model.beta;
    ll = N*model.outputDim*(log(model.beta) - log(2*pi)) - sum(sum(centred.*centred))* ...
         model.beta;
    ll = ll*0.5;
def model = linearOptimise(model, X, Y, varargin):

    """Optimise a linear model.
    
    Description:
    
    model = linearOptimise(model, X, Y) optimises a linear model using
     least squares.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
      X - the input data locations for the optimisation.
      Y - the target data locations for the optimisation.
        

    See also
    LINEARCREATE, MODELOPTIMISE


    Copyright (c) 2005, 2006, 2008 Neil D. Lawrence
    
    """
        
    N = size(X, 1);
    Xo = [X ones(N, 1)];
    W = pdinv(Xo'*Xo)*Xo'*Y;
    model.b = W(end, :);
    model.W = W(1:end-1, :);
    %model.b = mean(Y); %W(end, :);
    %model.W = pdinv(X'*X)*X'*(Y - repmat(model.b, size(Y, 1), 1));
    centred = Y - linearOut(model, X);
    model.beta = N./sum(centred.*centred, 1);

def options = linearOptions;

    """Options for learning a linear model.
    
    Description:
    def options = linearOptions;
%
    """
        
    options.activeFunc = 'linear';
def Y = linearOut(model, X):

    """Obtain the output of the linear model.
    
    Description:
    
    Y = linearOut(model, X) gives the output of a linear model.
     Returns:
      Y - the output.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
        

    See also
    MODELOUT, LINEARCREATE


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    numData = size(X, 1);
    Y = X*model.W + ones(numData, 1)*model.b;

def g = linearOutputGrad(model, X):

    """Evaluate derivatives of linear model outputs with respect to parameters.
    
    Description:
    def g = linearOutputGrad(model, X):
%
    """
        
    numData = size(X, 1);
    for i = 1:model.outputDim
      startZeros = zeros(numData, model.inputDim*(i - 1));
      finishZeros = zeros(numData, model.inputDim*(model.outputDim-i));
      startZeros2 = zeros(numData, (i - 1));
      finishZeros2 = zeros(numData, (model.outputDim-i));
      g(:, :, i) = [startZeros X finishZeros startZeros2 ones(numData, 1) finishZeros2];
    end
def g = linearOutputGradX(model, X):

    """Evaluate derivatives of linear model outputs with respect to inputs.
    
    Description:
    
    g = linearOutputGradX(model, X) returns the derivatives of the
     outputs of an LINEAR model with respect to the inputs to the
     model.
     Returns:
      g - the gradient of the output with respect to the inputs, in
     Arguments:
      model - the model for which the derivatives will be computed.
      X - the locations at which the derivatives will be computed.
        

    See also
    LINEAROUTPUTGRAD, MODELOUTPUTGRADX


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    g = repmat(shiftdim(model.W, -1), [size(X, 1) 1 1]);
def model = linearParamInit(model):

    """Initialise the parameters of an LINEAR model.
    
    Description:
    
    model = linearParamInit(model) sets the initial weight vectors and
     biases to small random values.
     Returns:
      model - the initialised model.
     Arguments:
      model - the input model to initialise.
        

    See also
    MODELPARAMINIT, LINEARCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    model.W = randn(model.inputDim, model.outputDim)/sqrt(model.inputDim + 1);
    model.b = randn(1, model.outputDim)/sqrt(model.inputDim + 1);
    model.beta = 1;
def model = lleCreate(inputDim, outputDim, Y, options):

    """Locally linear embedding model.
    
    Description:
    
    model = lleCreate(inputDimension, outputDim, Y, options) creates a
     structure for a locally linear embedding.
     Returns:
      model - model structure containing the neural network specified.
     Arguments:
      inputDimension - dimension of latent space.
      outputDim - dimension of data.
      Y - the data to be modelled in design matrix format (as many rows
       as there are data points).
      options - options structure as returned by lleCreate.
        

    See also
    LLEOPTIONS, MODELCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    
    model.type = 'lle';
    
    if size(Y, 2) ~= outputDim
      error(['Input matrix Y does not have dimension ' num2str(d)]);
    end
    model.k = options.numNeighbours;
    model.Y = Y;
    model.d = outputDim;
    model.q = options.latentDim;
    model.N = size(Y, 1);
def X = lleEmbed(Y, dims, neighbours):

    """Embed data set with LLE.
    
    Description:
    def X = lleEmbed(Y, dims, neighbours):
%
    """
        
    % Wrapper for Sam Roweis' LLE code.
    
    % Note LLE code uses the transpose of a design matrix.
    if nargin < 3
      neighbours = 7;
    end
    if any(any(isnan(Y)))
      error('Cannot initialise gplvm using LLE when missing data is present.')
    else
      X = lle(Y', neighbours, dims);
      X = X';
      % Rescale X so that variance is 1 and mean is zero.
      meanX = mean(X);
      X = X-ones(size(Y, 1), 1)*meanX;
      varX = var(X);
      X = X*diag(sqrt(1./varX));
    end

def model = lleOptimise(model, display, iters):

    """Optimise an LLE model.
    
    Description:
    
    model = lleOptimise(model) optimises an mixtures of Gaussians
     model via the expectation maximisation algorithm.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
        

    See also
    MMPCACREATE, MODELOPTIMISE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    model.indices = findNeighbours(model.Y, model.k);
    model.W = spalloc(model.N, model.N, model.N*model.k);
    for i = 1:model.N
      Ytemp = model.Y(model.indices(i, :), :);
      %  C = cov([model.Y(i, :); Ytemp]');
      Ytemp = Ytemp - repmat(model.Y(i, :), model.k, 1);
      C = Ytemp*Ytemp';
      if model.d<model.k
        %C = C + trace(C)*1e-3*eye(model.k+1);
        C = C + trace(C)*1e-3*eye(model.k);
      end
      %[fullU, v] = eig(C);
      %[void, ind] = min(diag(v));
      %what = fullU(:, ind);
      [U, jitter] = jitChol(C);
      y = U'\ones(model.k, 1);
      what = U\y;
      
      %what = (C\ones(model.k, 1));
      what = what/sum(what);
      model.W(i, model.indices(i, :)) = what';
      model.W(i, i) = -1;
      %    model.W(i, [i model.indices(i, :)]) = what';
    end
    
    options.disp = 0; 
    options.isreal = 1; 
    options.issym = 1; 
    %[m, v] = svds(model.W', model.q+1, 0);
    
    %model.P = speye(model.N) - model.W - model.W' +model.W'*model.W;
    model.P = model.W'*model.W;
    if isoctave
      warning('No eigs function in Octave');
      % Nasty hack for eigenvalue problem in Octave.
      [m, v] = eig(model.P);
      [v, order] = sort(diag(v));
      v = diag(v(1:model.q+1));
      m = m(:, order);
      m = m(:, 1:model.q+1);
    else
      [m, v] = eigs_r11(model.P, model.q+1, 'sm', options);
    end
    [void, ind] = min(diag(v));
    model.X = m(:, [1:(ind-1) (ind+1):end]);
def options = lleOptions(neighbours, latentDim):

    """Options for a density network.
    
    Description:
    
    options = lleOptions returns the default options for a locally
     linear embedding.
     Returns:
      options - default options structure for locally linear embedding.
        

    See also
    LLECREATE, MLPCREATE, RBFCREATE, KBRCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    if nargin < 2
      latentDim = 2
      if nargin < 1
        neighbours = 7;
      end
    end
    options.latentDim = latentDim;
    options.numNeighbours = neighbours;
    
    

def X = lmvuEmbed(Y,dims,k,nr_landmark):

    """Embed data set with landmark MVU
    
    Description:
    
    X = lmvuEmbed(Y, dims, k, nr_landmark) Embed data set with
     landmark version of Weinberg et al.'s maximum variance unfolding
     algorithm.
     Returns:
      X - embedding
     Arguments:
      Y - Data
      dims - Dimensionality of Embedding (default = 2)
      k - Number of Neighbours in Proximity Graph (default = 7)
      nr_landmark - Number of landmark Points
        

    See also
    PPCAEMBED, LLEEMBED, MVUEMBED


    Copyright (c) Neil D. Lawrence, 2007 Carl Henrik Ek
    
    """
        
    if(nargin<4)
      nr_landmark = 30;
      if(nargin<3)
        k = 7;
        if(nargin<2)
          dims = 2;
          if(nargin<1)
    	error('To Few Arguments');
          end
        end
      end
    end
      
    if(any(any(isnan(Y))))
      error('Cannot Initialise GPLVM using lmvu when missing data is present.');
    end
    
    X = lmvu(distance(Y'),nr_landmark,k);
    
    X = X(1:1:dims,:)';
    
    return
def lvmClassVisualise(call):

    """Callback function for visualising data in 2-D.
    
    Description:
    def lvmClassVisualise(call):
%
    """
        
    global visualiseInfo
    
    
    switch call
     case 'click'
      [x, y]  = localCheckPointPosition(visualiseInfo);  
      if ~isempty(x) 
        visualiseInfo.latentPos = [x, y];
      end
      visualiseInfo.clicked = ~visualiseInfo.clicked;
      if isfield(visualiseInfo.model, 'dynamics') & ~isempty(visualiseInfo.model.dynamics)
        if visualiseInfo.runDynamics
          visualiseInfo.dynamicsRunning = 1;
          fhandle = str2func([visualiseInfo.model.type 'DynamicsRun']);
          feval(fhandle);
          visualiseInfo.dynamicsRunning = 0;
        end
      else
        visualiseInfo.dynamicsRunning = 0;
      end
     case 'move'
      if visualiseInfo.clicked & ~visualiseInfo.runDynamics
        [x, y]  = localCheckPointPosition(visualiseInfo);  
        if ~isempty(x) 
          % This should be changed to a model specific visualisation.
          set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);
          fhandle = str2func([visualiseInfo.model.type 'PosteriorMeanVar']);
          [mu, varsigma] = fhandle(visualiseInfo.model, [x y]);
          if isfield(visualiseInfo.model, 'noise')
            Y = noiseOut(visualiseInfo.model.noise, mu, varsigma);
          else
            Y = mu;
          end
          visualiseInfo.visualiseModify(visualiseInfo.visHandle, ...
                                        Y, visualiseInfo.varargin{:});
          visualiseInfo.latentPos=[x, y];
        end
      end
     case 'toggleDynamics'
      visualiseInfo.runDynamics = ~visualiseInfo.runDynamics;
      set(visualiseInfo.dynamicsRadio, 'value', visualiseInfo.runDynamics);
    
     case 'dynamicsSliderChange'
      X = modelOut(visualiseInfo.model.dynamics, get(visualiseInfo.dynamicsSlider, 'value'));
      x = X(1);
      y = X(2);
      visualiseInfo.latentPos = [x, y];
      set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);
      fhandle = str2func([visualiseInfo.model.type 'PosteriorMeanVar']);
      [mu, varsigma] = fhandle(visualiseInfo.model, [x y]);
      if isfield(visualiseInfo.model, 'noise')
        Y = noiseOut(visualiseInfo.model.noise, mu, varsigma);
      else
        Y = mu;
      end
      visualiseInfo.visualiseModify(visualiseInfo.visHandle, ...
                                    Y, visualiseInfo.varargin{:});
    end
    
    
    
    
def point = localGetNormCursorPoint(figHandle):
    
    point = get(figHandle, 'currentPoint');
    figPos = get(figHandle, 'Position');
    % Normalise the point of the curstor
    point(1) = point(1)/figPos(3);
    point(2) = point(2)/figPos(4);
    
def [x, y] = localGetNormAxesPoint(point, axesHandle):
    
    position = get(axesHandle, 'Position');
    x = (point(1) - position(1))/position(3);
    y = (point(2) - position(2))/position(4);
    lim = get(axesHandle, 'XLim');
    x = x*(lim(2) - lim(1));
    x = x + lim(1);
    lim = get(axesHandle, 'YLim');
    y = y*(lim(2) - lim(1));
    y = y + lim(1);
    
    
def [x, y] = localCheckPointPosition(visualiseInfo):
    
    % Get the point of the cursor
    point = localGetNormCursorPoint(gcf);
    
    % get the position of the axes
    position = get(visualiseInfo.plotAxes, 'Position');
    
    
    % Check if the pointer is in the axes
    if point(1) > position(1) ...
          & point(1) < position(1) + position(3) ...
          & point(2) > position(2) ...
          & point(2) < position(2) + position(4);
      
      % Rescale the point according to the axes
      [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);
    
      % Find the nearest point
    else
      % Return nothing
      x = [];
      y = [];
    end

def lvmClassVisualisePath(call):

    """Latent variable model path drawing in latent space.
    
    Description:
    def lvmClassVisualisePath(call):
%
    """
        
    global visualiseInfo
    
    
    switch call
     case 'click'
      visualiseInfo.clicked = ~visualiseInfo.clicked;
      if visualiseInfo.clicked
        [x, y]  = localCheckPointPosition(visualiseInfo);  
        if ~isempty(x) 
          visualiseInfo.latentPos = [x, y];
          visualiseInfo.posTrail = [x, y];
          visualiseInfo.lastTime = clock;
          visualiseInfo.timeTrail = 0;
        end
      end
     case 'move'
      if visualiseInfo.clicked 
        [x, y]  = localCheckPointPosition(visualiseInfo);  
        if ~isempty(x) 
          visualiseInfo.latentPos=[x, y];
          % Draw connecting lines
          iplot(visualiseInfo.plotAxes, x, y);
          visualiseInfo.posTrail = [visualiseInfo.posTrail; [x, y]];
          timeNow = clock;
          visualiseInfo.timeTrail = [visualiseInfo.timeTrail; visualiseInfo.timeTrail(end)+etime(timeNow, ...
                                                            visualiseInfo.lastTime)];
          visualiseInfo.lastTime = timeNow;
        end
      end
    end
    
    
    
    
def point = localGetNormCursorPoint(figHandle):
    
    point = get(figHandle, 'currentPoint');
    figPos = get(figHandle, 'Position');
    % Normalise the point of the curstor
    point(1) = point(1)/figPos(3);
    point(2) = point(2)/figPos(4);
    
def [x, y] = localGetNormAxesPoint(point, axesHandle):
    
    position = get(axesHandle, 'Position');
    x = (point(1) - position(1))/position(3);
    y = (point(2) - position(2))/position(4);
    lim = get(axesHandle, 'XLim');
    x = x*(lim(2) - lim(1));
    x = x + lim(1);
    lim = get(axesHandle, 'YLim');
    y = y*(lim(2) - lim(1));
    y = y + lim(1);
    
    
def [x, y] = localCheckPointPosition(visualiseInfo):
    
    % Get the point of the cursor
    point = localGetNormCursorPoint(gcf);
    
    % get the position of the axes
    position = get(visualiseInfo.plotAxes, 'Position');
    
    
    % Check if the pointer is in the axes
    if point(1) > position(1) ...
          & point(1) < position(1) + position(3) ...
          & point(2) > position(2) ...
          & point(2) < position(2) + position(4);
      
      % Rescale the point according to the axes
      [x y] = localGetNormAxesPoint(point, visualiseInfo.plotAxes);
    
      % Find the nearest point
    else
      % Return nothing
      x = [];
      y = [];
    end

def [model, lbls] = lvmLoadResult(modelType, dataSet, number):

    """Load a previously saved result.
    
    Description:
    
    [model, lbls] = lvmLoadResult(dataSet, number) loads a previously
     saved LVM result.
     Returns:
      model - the saved model.
      lbls - labels of the data set (for visualisation purposes).
     Arguments:
      dataSet - the name of the data set to load.
      number - the number of the LVM run to load.
        

    See also
    LVMLOADDATA


    Copyright (c) 2003, 2004, 2005, 2006, 2008 Neil D. Lawrence
    
    """
        
    
    
    [Y, lbls] = lvmLoadData(dataSet);
    
    dataSet(1) = upper(dataSet(1));
    if ~isempty(modelType)
        modelType(1) = upper(modelType(1));
    end
    
    load(['dem' dataSet modelType num2str(number)])

def err = lvmNearestNeighbour(model, lbls):

    """Give the number of errors in latent space for 1 nearest neighbour.
    
    Description:
    
    lvmNearestNeighbour(model, lbls) computes the number errors for 1
     nearest neighbour in latent space.
     Arguments:
      model - the model for which the computation is required.
      lbls - the labels of the data.


    Copyright (c) 2004, 2006, 2008 Neil D. Lawrence
    
    """
        
    d = dist2(model.X, model.X);
    for i = 1:size(model.X, 1); 
      d(i, i) = inf; 
    end
    
    for i= 1:size(lbls, 1); 
      lbls2(i, :) =  find(lbls(i, :));
    end
    [void, ind] = min(d);
    err = size(model.X, 1) - sum(lbls2(ind) == lbls2);

def lvmPrintPlot(model, lbls, capName, experimentNo, colour):

    """Print latent space for learnt model.
    
    Description:
    
    lvmPrintPlot(model, lbls, capName, experimentNo) prints a latent
     space repsresentation for an LVM model.
     Arguments:
      model - the model to use for plotting the latent space.
      lbls - any lables that are available for plotting.
      capName - the name of the saved plots.
      experimentNo - the experiment number to assign to the files.
        

    See also
    LVMSCATTERPLOT


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if nargin < 5
      colour = 0;
    end
    if colour
      lvmScatterPlotColor(model, lbls);
    else
      lvmScatterPlot(model, lbls);
    end
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    fileName = ['dem' capName modelType num2str(experimentNo)];
    directory = ['../tex/diagrams'];
    printPlot(fileName, directory, '../html');
    
    figure
    clf
    ax = axes('position', [0.05 0.05 0.9 0.9]);
    hold on
    if ~isempty(lbls) && ~strcmp(lbls, 'connect')
      lvmTwoDPlot(model.X, lbls, getSymbols(size(lbls, 2)));
    else
      lvmTwoDPlot(model.X, lbls);
    end
    xLim = [min(model.X(:, 1)) max(model.X(:, 1))]*1.1;
    yLim = [min(model.X(:, 2)) max(model.X(:, 2))]*1.1;
    set(ax, 'xLim', xLim);
    set(ax, 'yLim', yLim);
    
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    printPlot([fileName 'NoGray'], directory, '../html');
    %print('-depsc', ['../tex/diagrams/' fileName 'NoGray'])
    %print('-deps', ['../tex/diagrams/' fileName 'NoGrayNoColour'])

def lvmResultsDynamic(modelType, dataSet, number, dataType, varargin):

    """Load a results file and visualise them.
    
    Description:
    
    lvmResultsDynamic(modelType, dataSet, number, dataType, ...) loads
     results of a latent variable model and visualises them.
     Arguments:
      modelType - the type of model ran on the data set.
      dataSet - the name of the data set to load.
      number - the number of the run used.
      dataType - the type of data to visualise.
      ... - additional arguments to be passed to the lvmVisualise
       command.
        

    See also
    LVMLOADRESULT, LVMVISUALISE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    [model, lbls] = lvmLoadResult(modelType, dataSet, number);
    
    % Visualise the results
    switch size(model.X, 2) 
     case 2
      lvmVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                     varargin{:});
      
     otherwise 
      error('No visualisation code for data of this latent dimension.');
    end
def [ax, data] = lvmScatterPlot(model, YLbls, ax):

    """2-D scatter plot of the latent points.
    
    Description:
    
    ax = lvmScatterPlot(model) produces a visualisation of the latent
     space with the given model.
     Returns:
      ax - the axes handle where the scatter plot was placed.
     Arguments:
      model - the model for which the scatter plot is being produced.
        DESC produces a visualisation of the latent space for the given model,
        using the provided labels to distinguish the latent points.
        ARG model : the model for which the scatter plot is being produced.
        ARG lbls : labels for each data point so that they may be given different
        symbols. Useful when each data point is associated with a different
        class.
        RETURN ax : the axes handle where the scatter plot was placed.
        
        DESC produces a visualisation of the latent space for the given model,
        using the provided labels to distinguish the latent points.
        ARG model : the model for which the scatter plot is being produced.
        ARG lbls : labels for each data point so that they may be given different
        symbols. Useful when each data point is associated with a different
        class.
        ARG ax : the axes where the plot is to be placed.
        RETURN ax : the axes handle where the scatter plot was placed.
        
        

    See also
    LVMVISUALISE, LVMTWODPLOT, LVMSCATTERPLOTCOLOR


    Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
    
    """
        
    if nargin<3
      ax = [];
      if nargin < 2
        YLbls = [];
      end
    end
    
    if isempty(YLbls)
      symbol = [];
    else
      symbol = getSymbols(size(YLbls,2));
    end
    
    x1Min = min(model.X(:, 1));
    x1Max = max(model.X(:, 1));
    x1Span = x1Max - x1Min;
    x1Min = x1Min - 0.05*x1Span;
    x1Max = x1Max + 0.05*x1Span;
    x1 = linspace(x1Min, x1Max, 150);
    
    x2Min = min(model.X(:, 2));
    x2Max = max(model.X(:, 2));
    x2Span = x2Max - x2Min;
    x2Min = x2Min - 0.05*x2Span;
    x2Max = x2Max + 0.05*x2Span;
    x2 = linspace(x2Min, x2Max, 150);
    
    funcStr = [model.type 'PosteriorMeanVar'];
    if exist(funcStr)==2
      [X1, X2] = meshgrid(x1, x2);
      XTest = [X1(:), X2(:)];
      fhandle = str2func(funcStr);
      if str2num(version('-release'))>13
        [mu, varsigma] = fhandle(model, XTest);
      else 
        [mu, varsigma] = feval(fhandle, model, XTest);
      end
      
      d = size(mu, 2);
      if size(varsigma, 2) == 1
        dataMaxProb = -0.5*d*log(varsigma);
      else
        dataMaxProb = -.5*sum(log(varsigma), 2);
      end
      
      if isempty(ax)
        figure(1)
        clf
        % Create the plot for the data
        ax = axes('position', [0.05 0.05 0.9 0.9]);
      else
        axes(ax);
      end
      hold on
    
      C = reshape(dataMaxProb, size(X1));
      
      % Rescale it
      C = C - min(min(C));
      if max(max(C))~=0
        C = C/max(max(C));
        C = round(C*63);
        image(x1, x2, C);
      end
      
      %[c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
      % shading flat
      colormap gray;
      %colorbar
    end
    data = lvmTwoDPlot(model.X, YLbls, symbol);
    switch model.type
     case 'dnet'
      plot(model.X_u(:, 1), model.X_u(:, 2), 'g.')
    end
    xLim = [min(x1) max(x1)];
    yLim = [min(x2) max(x2)];
    set(ax, 'xLim', xLim);
    set(ax, 'yLim', yLim);
    
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    

def [ax, data] = lvmScatterPlotColor(model, shade, ax):

    """2-D scatter plot of the latent points with color.
    
    Description:
    
    ax = lvmScatterPlotColor(model, shade) produces a visualisation of
     the latent space with the given model using color shadings given,
     specifically written for the 'swiss roll data'.
     Returns:
      ax - the axes handle where the scatter plot was placed.
     Arguments:
      model - the model for which the scatter plot is being produced.
      shade - color indicator for each data point so that they may be
       given different colours.
        DESC produces a visualisation of the latent space for the given model,
        using the provided labels to distinguish the latent points.
        ARG model : the model for which the scatter plot is being produced.
        ARG shade : colour indicator for each data point so that they may be given
        different colours.
        ARG ax : the axes where the plot is to be placed.
        RETURN ax : the axes handle where the scatter plot was placed.
        
        

    See also
    FGPLVMVISUALISE, LVMTWODPLOT, LVMSCATTERPLOT


    Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence
    
    """
        
    if nargin < 3
      ax = [];
    end
    
    shade = shade - min(shade)+eps;
    shade = shade/max(shade);
    shade = ceil(shade*64);
    
    x1 = linspace(min(model.X(:, 1))*1.1, max(model.X(:, 1))*1.1, 30);
    x2 = linspace(min(model.X(:, 2))*1.1, max(model.X(:, 2))*1.1, 30);
    
    funcStr = [model.type 'PosteriorMeanVar'];
    if exist(funcStr)==2
      [X1, X2] = meshgrid(x1, x2);
      XTest = [X1(:), X2(:)];
      
      fhandle = str2func([model.type 'PosteriorMeanVar']);
      if str2num(version('-release'))>13
        [mu, varsigma] = fhandle(model, XTest);
      else 
        [mu, varsigma] = feval(fhandle, model, XTest);
      end
      
      d = size(mu, 2);
      if size(varsigma, 2) == 1
        dataMaxProb = -0.5*d*log(varsigma);
      else
        dataMaxProb = -.5*sum(log(varsigma), 2);
      end
      
      if isempty(ax)
        figure(1)
        clf
        % Create the plot for the data
        ax = axes('position', [0.05 0.05 0.9 0.9]);
      else
        axes(ax);
      end
      hold on
    
      C = reshape(dataMaxProb, size(X1));
      
      % Rescale it
      C = C - min(min(C));
      if max(max(C))~=0
        C = C/max(max(C));
        C = round(C*63);
        image(x1, x2, C);
      end
      
      [c, h] = contourf(X1, X2, log10(reshape(1./varsigma(:, 1), size(X1))), 128); 
      colormap gray;
      
      
      figure(1)
      clf
      shading flat
      gr = colormap('gray');
      set(h, 'CDataMapping', 'direct')
      for i = 1:length(h)
        set(h(i), 'cdata', i);
      end
      %colorbar
    end
    jt = colormap('jet');
    data = lvmTwoDPlot(model.X);
    
    for i=1:length(data)
      set(data(i), 'color', jt(shade(i), :));
    end
    xLim = [min(x1) max(x1)];
    yLim = [min(x2) max(x2)];
    set(ax, 'xLim', xLim);
    set(ax, 'yLim', yLim);
    
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    

def [ax, data] = lvmScatterPlotNeighbours(model, YLbls, ax):

    """2-D scatter plot of the latent points with neighbourhood.
    
    Description:
    
    ax = lvmScatterPlotNeighbours(model) produces a visualisation of
     the latent space with the given model.
     Returns:
      ax - the axes handle where the scatter plot was placed.
     Arguments:
      model - the model for which the scatter plot is being produced.
        DESC produces a visualisation of the latent space for the given model,
        using the provided labels to distinguish the latent points.
        ARG model : the model for which the scatter plot is being produced.
        ARG lbls : labels for each data point so that they may be given different
        symbols. Useful when each data point is associated with a different
        class.
        RETURN ax : the axes handle where the scatter plot was placed.
        
        DESC produces a visualisation of the latent space for the given model,
        using the provided labels to distinguish the latent points.
        ARG model : the model for which the scatter plot is being produced.
        ARG lbls : labels for each data point so that they may be given different
        symbols. Useful when each data point is associated with a different
        class.
        ARG ax : the axes where the plot is to be placed.
        RETURN ax : the axes handle where the scatter plot was placed.
        
        

    See also
    LVMSCATTERPLOT, LVMTWODPLOT, LVMSCATTERPLOTCOLOR


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    if nargin<3
      ax = [];
      if nargin < 2
        YLbls = [];
      end
    end
    
    if isempty(YLbls)
      symbol = [];
    else
      symbol = getSymbols(size(YLbls,2));
    end
    
    [ax, data] = lvmScatterPlot(model, YLbls, ax);
    
    
    for i = 1:model.N
      for j = 1:length(model.indices(i, :))
        line([model.X(i, 1)  ...
             model.X(model.indices(i, j), 1)], ...
             [model.X(i, 2) ...
             model.X(model.indices(i, j), 2)]);
      end
    end

def returnVal = lvmTwoDPlot(X, lbl, symbol):

    """Helper function for plotting the labels in 2-D.
    
    Description:
    
    lvmTwoDPlot(X, lbl, symbol) helper function for plotting an
     embedding in 2-D with symbols.
     Arguments:
      X - the data to plot.
      lbl - the labels of the data point.
      symbol - the symbols to use for the different labels.
        

    See also
    LVMSCATTERPLOT, LVMVISUALISE


    Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence
    
    """
        if nargin < 2
      lbl = [];
    end
    if(strcmp(lbl, 'connect'))
      connect = true;
      lbl = [];
    else
      connect = false;
    end
    
    if nargin < 3
      if isempty(lbl)
        symbol = getSymbols(1);
      else
        symbol = getSymbols(size(lbl,2));
      end
    end
    axisHand = gca;
    returnVal = [];
    nextPlot = get(axisHand, 'nextplot');
    for i = 1:size(X, 1)
      if i == 2
        set(axisHand, 'nextplot', 'add');
      end
      if ~isempty(lbl)
        labelNo = find(lbl(i, :));
      else
        labelNo = 1;
      end
      try 
        returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
        if connect
          if i>1
            line([X(i-1, 1) X(i, 1)], [X(i-1, 2) X(i, 2)]);
          end
        end
      catch
        if strcmp(lasterr, 'Index exceeds matrix dimensions.')
          error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
        end
      end
    end
    set(axisHand, 'nextplot', nextPlot);
    set(returnVal, 'markersize', 10);
    set(returnVal, 'linewidth', 2);

def lvmVisualise(model, YLbls, ...
			visualiseFunction, visualiseModify, varargin):

    """Visualise the manifold.
    
    Description:
    
    lvmVisualise(model, Ylbls, visualiseFunction, visualiseModify,
     ...) visualises a two dimensional manifold in data space using
     commands passed as argument.
     Arguments:
      model - the model to visualise (of type lvm).
      Ylbls - any labels for the training data to improve the
       visualisation.
      visualiseFunction - the function that draws the visualisation (in
       data space) when the graphs are first drawn.
      visualiseModify - the function that modifies the visualisation as
       you move around the latent space.
      ... - various additional arguments to be passed to the
       visualisation commands.
        

    See also
    LVMSCATTERPLOT, LVMRESULTSDYNAMIC


    Copyright (c) 2003, 2004, 2005, 2006, 2008 Neil D. Lawrence
    
    """
        
    global visualiseInfo
    
    figure(1)
    clf
    visualiseInfo.plotAxes = lvmScatterPlot(model, YLbls);
    set(get(visualiseInfo.plotAxes, 'title'), 'string', 'X', 'fontsize', 30);
    set(visualiseInfo.plotAxes, 'position', [0.05 0.05 0.9 0.8]);
    
    visualiseInfo.latentHandle = line(0, 0, 'markersize', 20, 'color', ...
                                      [0 0 0], 'marker', '.', 'visible', ...
                                      'on', 'erasemode', 'xor');
    
    % Set up the X limits and Y limits of the main plot
    xLim = [min(model.X(:, 1)) max(model.X(:, 1))];
    xSpan = xLim(2) - xLim(1);
    xLim(1) = xLim(1) - 0.05*xSpan;
    xLim(2) = xLim(2) + 0.05*xSpan;
    xSpan = xLim(2) - xLim(1);
    
    yLim = [min(model.X(:, 2)) max(model.X(:, 2))];
    ySpan = yLim(2) - yLim(1);
    yLim(1) = yLim(1) - 0.05*ySpan;
    yLim(2) = yLim(2) + 0.05*ySpan;
    ySpan = yLim(2) - yLim(1);
    
    set(visualiseInfo.plotAxes, 'XLim', xLim)
    set(visualiseInfo.plotAxes, 'YLim', yLim)
    
    visualiseInfo.clicked = 0;
    
    visualiseInfo.digitAxes = [];
    visualiseInfo.digitIndex = [];
    
    visualiseInfo.dynamicsRadio = ...
        uicontrol('Style', 'radiobutton', ...
                  'String', 'Run Dynamics', ...
                  'units', 'normalized', ...
                  'position', [0 0 0.2 0.05], ...
                  'Callback', 'lvmClassVisualise(''toggleDynamics'')', ...
                  'value', 0);
    
    visualiseInfo.dynamicsSlider = ...
        uicontrol('Style', 'slider', ...
                  'String', 'Time', ...
                  'sliderStep', [0.01, 0.1], ...
                  'units', 'normalized', ...
                  'position', [0 0.95 1 0.05], ...
                  'callback', 'lvmClassVisualise(''dynamicsSliderChange'')');
    
    if ~isfield(model, 'dynamics') | isempty(model.dynamics)
      set(visualiseInfo.dynamicsRadio, 'visible', 'off');
      set(visualiseInfo.dynamicsSlider, 'visible', 'off');
    else
      if ~isfield(model.dynamics, 'dynamicsType') 
        set(visualiseInfo.dynamicsRadio, 'visible', 'on');
        set(visualiseInfo.dynamicsSlider, 'visible', 'off');
      else
        switch model.dynamics.dynamicsType
         case 'regressive'
          set(visualiseInfo.dynamicsRadio, 'visible', 'off');
          set(visualiseInfo.dynamicsSlider, 'visible', 'on');
          set(visualiseInfo.dynamicsSlider, 'min', min(model.dynamics.X), ...
                            'max', max(model.dynamics.X), ...
                            'value', model.dynamics.X(1))
         case 'auto-regressive'
          set(visualiseInfo.dynamicsRadio, 'visible', 'on');
          set(visualiseInfo.dynamicsSlider, 'visible', 'off');
        end
      end
    end
    visualiseInfo.runDynamics = false;
    % Set the callback function
    set(gcf, 'WindowButtonMotionFcn', 'lvmClassVisualise(''move'')')
    set(gcf, 'WindowButtonDownFcn', 'lvmClassVisualise(''click'')')
    
    figure(2)
    clf
    
    if length(visualiseFunction)>4 & strcmp(visualiseFunction(1:5), 'image') & length(varargin)>0
      set(gcf, 'menubar', 'none')
      xPixels = 115;
      yPixels = 115;
      set(gcf, 'position', [232 572 xPixels yPixels/varargin{1}(1)*varargin{1}(2)])
      visualiseInfo.visualiseAxes = subplot(1, 1, 1);
      xWidth = varargin{1}(1)/xPixels;
      yHeight = varargin{1}(2)/yPixels;
      set(visualiseInfo.visualiseAxes, 'position', [0.5-xWidth/2 0.5-yHeight/2 xWidth yHeight])
    else
      visualiseInfo.visualiseAxes =subplot(1, 1, 1);
    end
    visData = zeros(1,model.d);
    if(length(visualiseFunction)>4 & strcmp(visualiseFunction(1:5), 'image'))
      visData(1) = min(min(model.y));
      visData(end) = max(max(model.y));
    else
      [void, indMax]= max(sum((model.y.*model.y), 2));
      visData = model.y(indMax, :);
    end
    visualiseInfo.visualiseFunction = str2func(visualiseFunction);
    visHandle = visualiseInfo.visualiseFunction(visData, varargin{:});
    set(visHandle, 'erasemode', 'xor')
    colormap gray
    
    % Pass the data to visualiseInfo
    visualiseInfo.model = model;
    visualiseInfo.varargin = varargin;
    visualiseInfo.visualiseModify = str2func(visualiseModify);
    visualiseInfo.visHandle = visHandle;
    
    set(get(visualiseInfo.visualiseAxes, 'title'), 'string', 'Y', 'fontsize', ...
                      30);
    set(visualiseInfo.visualiseAxes, 'position', [0.05 0.05 0.9 0.8]);
    
    hold off
    
    
    
    

def model = mapmodelReadFromFID(FID, varargin):

    """Load from a FID produced by C++ code.
    
    Description:
    
    model = mapmodelReadFromFID(FID) loads in from a file stream the
     data format produced by C++ code.
     Returns:
      model - the model loaded in from the file.
     Arguments:
      FID - the file ID from where the data is loaded.
        

    See also
    MODELREADFROMFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    modelType = readStringFromFID(FID, 'type');
    feval = str2func([modelType 'ReadFromFID']);
    model = feval(FID, varargin{:});

def model = mappingOptimise(model, X, Y, varargin):

    """Optimise the given model.
    
    Description:
    def model = mappingOptimise(model, X, Y, varargin):
%
    """
        
    fhandle = str2func([model.type 'Optimise']);
    model = fhandle(model, X, Y, varargin{:});
def X = matrixReadFromFID(FID, varargin):

    """Read a matrix from an FID.
    
    Description:
    
    X = matrixReadFromFID(FID) reads a matrix from an FID.
     Returns:
      X - the returned matrix read from the file.
     Arguments:
      FID - the file ID to read the matrix from.
        

    See also
    MODELREADFROMFID, DOUBLEMATRIXREADFROMFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    modelType = readStringFromFID(FID, 'type');
    feval = str2func([modelType 'ReadFromFID']);
    X = feval(FID, varargin{:});

def model = mlpCreate(inputDim, outputDim, options):

    """Multi-layer peceptron model.
    
    Description:
    
    model = mlpCreate(inputDimension, outputDim, options) creates a
     structure for a multi-layer perceptron. For models with a single
     hidden layer it is a wrapper structure for NETLAB's multi-layer
     perceptron model.
     Returns:
      model - model structure containing the neural network specified.
     Arguments:
      inputDimension - dimension of input data.
      outputDim - dimension of target data.
      options - options structure. The structure contains the type of
       output 'activation function', the number of hidden units and the
       optimiser to be used. A set of default options are given by the
       file mlpOptions.
        

    See also
    MLPOPTIONS, MLP


    Copyright (c) 2005, 2006, 2007 Neil D. Lawrence
    
    """
        
    
    if length(options.hiddenDim) == 1
      % Can use NETLAB implementation.
      model = mlp(inputDim, options.hiddenDim, outputDim, ...
                  options.activeFunc);
      model.numParams = model.nwts;
      model.hiddenDim = options.hiddenDim;
      model.inputDim = inputDim;
      model.outputDim = outputDim;
    else
      error('Multiple hidden layer mlp not yet implemented')
      model.type = 'mlp';
      model.hiddenDim = options.hiddenDim;
      model.inputDim = inputDim;
      model.outputDim = outputDim;
      model.numParams = (model.inputDim+1)*model.hiddenDim(1);
      for i = 1:length(model.hiddenDim)-1
        model.numParams = model.numParams + (model.hiddenDim(i)+1)*model.hiddenDim(i+1);
      end
      model.numParams = model.numParams + (model.hiddenDim(end)+1)*model.outputDim;
    
      activationFunctions = {'linear', 'logistic', 'softmax'};
    
      if sum(strcmp(options.activeFunc, activationFunctions)) == 0
        error('Undefined output function.');
      else
        model.outfn = options.activeFunc;
      end
      model = mlpParamInit(model);
    end
    
    model.optimiser = options.optimiser;

def mlpDisplay(model, spacing):

    """Display the multi-layer perceptron model.
    
    Description:
    
    mlpDisplay(model, spacing) displaces the contents of a multi-layer
     perceptron model.
     Arguments:
      model - the model to be displayed.
      spacing - optional spacing to place before model display.
        

    See also
    MLPCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Multi-layer perceptron model:\n')
    fprintf(spacing);
    fprintf('  Input units: %d\n', model.inputDim);
    fprintf(spacing);
    fprintf('  Output units: %d\n', model.outputDim);
    if length(model.hiddenDim)==1
      fprintf(spacing);
      fprintf('  Hidden units: %d\n', model.hiddenDim);
    else
      fprintf(spacing);
      fprintf('  Hidden layers: %d\n', length(model.hiddenDim));
      for i = 1:length(model.hiddenDim)
        fprintf(spacing);
        fprintf('    Layer 1: %d nodes\n', model.hiddenDim(i));
      end
    end
    fprintf(spacing);
    fprintf('  Number of parameters: %d\n', model.numParams);
    fprintf(spacing);
    fprintf(['  Output function: ' model.outfn '\n']);
    

def model = mlpExpandParam(model, params):

    """Update mlp model with new vector of parameters.
    
    Description:
    
    model = mlpExpandParam(model, params) takes a vector of MLP
     weights and places them in their respective positions in the MLP
     model. For single hidden layer neural networks the function is a
     wrapper for the mlpunpak command.
     Returns:
      model - the model with the weights distributed in the correct
       places.
     Arguments:
      model - the model in which the weights are to be placed.
      params - a vector of the weights to be placed in the model.
        

    See also
    MLPUNPAK, MLPCREATE, MLPEXTRACTPARAM


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      model = mlpunpak(model, params);
    else
      startVal = 1;
      endVal = model.inputDim*model.hiddenDim(1);
      model.w{1} = reshape(params(startVal:endVal, model.inputDim, ...
                                model.hiddenDim(1)));
      startVal = endVal + 1;
      endVal = endVal + model.hiddenDim(1);
      model.b{1} = params(startVal:endVal);
      for i = 2:length(model.hiddenDim)
        startVal = endVal + 1;
        endVal = endVal + model.hiddenDim(i-1)*model.hiddenDim(i);
        model.w{i} = reshape(params(startVal:endVal), model.hiddenDim(i-1), ...
                                    model.hiddenDim(i));
        startVal = endVal + 1;
        endVal = endVal + model.hiddenDim(i);
        model.b{i} = params(startVal:endVal);
      end
      i = length(model.hiddenDim);
      startVal = endVal + 1;
      endVal = endVal + model.hiddenDim(i)*model.outputDim;
      model.w{i+1} = resphape(params(startVal:endVal), model.hiddenDim(i), ...
                              model.outputDim);
      startVal = endVal + 1;
      endVal = endVal + model.outputDim;
      model.b{i+1} = params(startVal:endVal);
    end

def [params, names] = mlpExtractParam(model):

    """Extract weights and biases from an MLP.
    
    Description:
    
    [params, names] = mlpExtractParam(model) returns a vector of all
     the weights and biases from a multi-layer perceptron model. For
     single hidden layer models the function is a wrapper for the
     mlppak command.
     Returns:
      params - vector of all the weights and biases returned by the
       model. The structure is governed by mlppak.
      names - optional additional returned cell array of the names of
       the parameters.
     Arguments:
      model - the model from which we wish to extract the weights and
       biases.
        

    See also
    MLPPAK, MLPCREATE, MLPEXPANDPARAM, MODELEXTRACTPARAM


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      params = mlppak(model);
      if nargout > 1
        counter = 0;
        for j = 1:size(model.w1, 2)
          for i = 1:size(model.w1, 1)
            counter = counter + 1;
            names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
          end
        end
        for j = 1:size(model.b1, 2)
          counter = counter + 1;
          names{counter} = ['Hidden node bias ' num2str(j)];
        end
        for j = 1:size(model.w2, 2)
          for i = 1:size(model.w2, 1)
            counter = counter + 1;
            names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
          end
        end
        for j = 1:size(model.b2, 2)
          counter = counter + 1;
          names{counter} = ['Output node bias ' num2str(j)];
        end
      end
    else
      params = zeros(1, model.numParams);
      startVal = 1;
      endVal = model.inputDim*model.hiddenDim(1);
      params(startVal:endVal) = model.w{1}(:)';
      startVal = endVal + 1;
      endVal = endVal + model.hiddenDim(1);
      params(startVal:endVal) = model.b{1};
      for i = 2:length(model.hiddenDim)
        startVal = endVal + 1;
        endVal = endVal + model.hiddenDim(i-1)*model.hiddenDim(i);
        params(startVal:endVal) = model.w{i}(:)';
        startVal = endVal + 1;
        endVal = endVal + model.hiddenDim(i);
        params(startVal:endVal) = model.b{i};
      end
      i = length(model.hiddenDim);
      startVal = endVal + 1;
      endVal = endVal + model.hiddenDim(i)*model.outputDim;
      params(startVal:endVal) = model.w{i+1}(:)';
      startVal = endVal + 1;
      endVal = endVal + model.outputDim;
      params(startVal:endVal) = model.b{i+1};
      if nargout > 1
        counter = 0;
        for j = 1:size(model.w{1}, 2)
          for i = 1:size(model.w{1}, 1)
            counter = counter + 1;
            names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
          end
        end
        for j = 1:size(model.b{1}, 2)
          counter = counter + 1;
          names{counter} = ['Hidden node bias ' num2str(j)];
        end
        for k = 2:length(model.hiddenDim)
          for j = 1:size(model.w{k}, 2)
            for i = 1:size(model.w{k}, 1)
              counter = counter + 1;
              names{counter} = ['Hidden weight layer ' num2str(k-1) '-' ...
                                num2str(k) ', node '  num2str(i) '-' ...
                                num2str(j)];
            end
          end
        end
        for j = 1:size(model.w{end}, 2)
          for i = 1:size(model.w{end}, 1)
            counter = counter + 1;
            names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
          end
        end
        for j = 1:size(model.b2, 2)
          counter = counter + 1;
          names{counter} = ['Output node bias ' num2str(j)];
        end
      end
    end
def g = mlpLogLikeGradients(model):

    """Multi-layer perceptron gradients.
    
    Description:
    
    g = mlpLogLikeGradients(model) computes the gradients of the log
     likelihood of a multi-layer perceptron with respect to the
     parameters. For networks with one hidden layer this is done by
     wrapping the mlpgrad command.
     Returns:
      g - the gradients of the model log likelihood.
     Arguments:
      model - the model structure for computing the log likelihood.
        

    See also
    MODELLOGLIKEIHOOD, MLPGRAD


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      g = -mlpgrad(model, model.X, model.y);
    else
      error('Not yet implemented.')
    end
def g = mlpLogLikeHessian(model):

    """Multi-layer perceptron Hessian.
    
    Description:
    
    g = mlpLogLikeHessian(model) computes the Hessian of the log
     likelihood of a multi-layer perceptron with respect to the
     parameters. For networks with a single hidden layer this is done
     by wrapping the mlpgrad command.
     Returns:
      g - the Hessian of the model log likelihood.
     Arguments:
      model - the model structure for computing the log likelihood.
        

    See also
    MODELLOGLIKEIHOOD, MLPGRAD


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim)==1
      g = -mlphess(model, model.X, model.y);
    else
      error('Hessian not yet available for this model.')
    end

def ll = mlpLogLikelihood(model):

    """Multi-layer perceptron log likelihood.
    
    Description:
    
    ll = mlpLogLikelihood(model) computes the log likelihood of a
     multi-layer perceptron model. For single hidden layer models this
     is done by wrapping the mlperr command.
     Returns:
      ll - the model log likelihood.
     Arguments:
      model - the model structure for computing the log likelihood.
        

    See also
    MODELLOGLIKEIHOOD, MLPERR


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    
    if length(model.hiddenDim) == 1
      ll = -mlperr(model, model.X, model.y);
    else
      Y = mlpOut(model, model.X);
      ll = -0.5*sum(sum((model.Y - Y).^2));
    end
    
    ll = ll - size(model.X, 1)/2*log(2*pi);
def model = mlpOptimise(model, X, Y, display, iters):

    """Optimise MLP for given inputs and outputs.
    
    Description:
    
    model = mlpOptimise(model, X, Y) optimises a MLP using a nonlinear
     optimiser. squares fit.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
      X - the input data locations for the optimisation.
      Y - the target data locations for the optimisation.
        

    See also
    MLPCREATE, MODELOPTIMISE


    Copyright (c) 2005, 2006, 2007 Neil D. Lawrence
    
    """
        
    
    if nargin < 4
      display = 1;
      if nargin < 5
        iters = 500;
      end
    end
    
    options = optOptions;
    options(14) = iters;
    options(1) = display;
    model = netopt(model, options, X, Y, model.optimiser);  
def options = mlpOptions(numHidden):

    """Options for the multi-layered perceptron.
    
    Description:
    
    options = mlpOptions returns the default options for a multi-layer
     perceptron.
     Returns:
      options - default options structure for Multi-layer peceptron.
    
    options = mlpOptions(numHidden)
     Returns:
      options - default options structure for Multi-layer peceptron with
       the specified number of hidden units.
     Arguments:
      numHidden - number of hidden units.
        

    See also
    MLPCREATE, MLP


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if nargin < 1
      numHidden = 20;
    end
    options.hiddenDim = numHidden;
    options.activeFunc = 'linear';
    options.optimiser = optimiDefaultOptimiser;
def [Y, G, A] = mlpOut(model, X):

    """Output of an MLP model.
    
    Description:
    
    Y = mlpOut(model, X) gives the output of a multi-layer perceptron
     model, for single hidden layer models the function is a wrapper
     for mlpfwd.
     Returns:
      Y - the output.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
    
    [Y, G] = mlpOut(model, X) gives the output of a multi-layer
     perceptron model.
     Returns:
      Y - the output.
      G - the hidden layer activations.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
        

    See also
    MLPFWD, MLP, MODELOUT


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim)==1
      if nargout > 1
        if nargout > 2
          [Y, G, A] = mlpfwd(model, X);
        else
          [Y, G] = mlpfwd(model, X);
        end
      else
        Y = mlpfwd(model, X);
      end
    else
      ndata = size(x, 1);
      G{1} = tanh(X*model.w{1} + repmat(model.b{1}, numData, 1));
      A{1} = G{1}*model.w{2} + repmat(model.b{2}, numData, 1);
      for i = 2:length(model.numHidden)
        G{i} = tanh(A{i-1});
        A{i} = G{i}*model.w{i+1} + repmat(model.b{i+1}, numData, 1);
      end
      switch model.outfn
        
       case 'linear' 
        y = A{end};
        
       otherwise 
        error('Output function not implemented in multiple hidden layer model.')
        
      end
    end    
    
      
def g = mlpOutputGrad(model, X):

    """Evaluate derivatives of mlp model outputs with respect to parameters.
    
    Description:
    
    g = mlpOutputGrad(model, X) evaluates the derivates of a
     multi-layer perceptron's outputs with respect to the parameters of
     the multi-layer perceptron. Currently it simply wraps the NETLAB
     mlpderiv function.
     Returns:
      g - the gradient of the outputs of the multi-layer perceptron with
       respect to each of the parameters. The size of the matrix is
       number of data x number of parameters x number of outputs of the
       model.
     Arguments:
      model - the model for which the derivatives are to be computed.
      X - the input data locations where the gradients are to be
       computed.
        

    See also
    MLPCREATE, MLPDERIV


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      g = mlpderiv(model, X);
    else
      
      [Y, G, A] = mlpOut(model, X);
      gw = cell(1, length(model.w));
      for i = 1:length(gw)
        gw{i} = zeros(size(model.w{i}));
      end   
      for i = 1:length(G);
        WdG{i} = (1-G{i}.*G{i})*w{i};
      end
      for k = 1:model.outputDim
        gw{end}(:, k) = Z{end}(:, k);
        gb{end} = 1;
        for i = length(model.w)-1:-1:1
          %gw{i} = 
          error('Not yet implemented');
        end
      end
    
    % Evaluate second-layer gradients.
    gw2 = z'*deltas;
    gb2 = sum(deltas, 1);
    
    % Now do the backpropagation.
    delhid = deltas*net.w2';
    delhid = delhid.*(1.0 - z.*z);
    
    % Finally, evaluate the first-layer gradients.
    gw1 = x'*delhid;
    gb1 = sum(delhid, 1);
    
    g = [gw1(:)', gb1, gw2(:)', gb2];
    end

def g = mlpOutputGradX(model, X):

    """Evaluate derivatives of mlp model outputs with respect to inputs.
    
    Description:
    
    g = mlpOutputGradX(model, X) returns the derivatives of the
     outputs of an MLP model with respect to the inputs to the model.
     Returns:
      g - the gradient of the output with respect to the inputs.
     Arguments:
      model - the model for which the derivatives will be computed.
      X - the locations at which the derivatives will be computed.
        

    See also
    MLPOUTPUTGRAD, MODELOUTPUTGRADX


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      [Y, Z] = mlpOut(model, X);
      gprime = 1-Z.*Z;
      g = zeros(size(X, 1), model.inputDim, model.outputDim);
      for i = 1:size(X, 1)
        for j = 1:size(X, 2)
          g(i, j, :) = shiftdim((gprime(i, :).*model.w1(j, :))*model.w2, -1);
        end
      end
    else
      error('Not yet implemented for more than one hidden layer')
    end

def model = mlpParamInit(model):

    """Initialise the parameters of an MLP model.
    
    Description:
    
    model = mlpParamInit(model) sets the initial weight vectors and
     biases to small random values.
     Returns:
      model - the initialised model.
     Arguments:
      model - the input model to initialise.
        

    See also
    MODELPARAMINIT, MLPCREATE


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    if length(model.hiddenDim) == 1
      model.w1 = randn(model.inputDim, model.nhidden)/sqrt(model.inputDim + 1);
      model.b1 = randn(1, model.nhidden)/sqrt(model.inputDim + 1);
      model.w2 = randn(model.nhidden, model.outputDim)/sqrt(model.nhidden + 1);
      model.b2 = randn(1, model.outputDim)/sqrt(model.nhidden + 1);
    else
      model.w{1} = randn(model.inputDim, model.hiddenDim(1))/sqrt(model.inputDim + 1);
      model.b{1} = randn(1, model.hiddenDim(1))/sqrt(model.inputDim + 1);
      for i = 2:length(model.hiddenDim)
        model.b{i} = randn(1, model.hiddenDim(i))/sqrt(model.hiddenDim(i-1) + 1);
        model.w{i} = randn(model.hiddenDim(i-1), model.hiddenDim(i))/sqrt(model.hiddenDim(i-1) + 1);
      end
      i = length(model.hiddenDim);
      model.b{i+1} = randn(1, model.outputDim)/sqrt(model.hiddenDim(i) + 1);
      model.w{i+1} = randn(model.hiddenDim(i), model.outputDim)/sqrt(model.hiddenDim(i) + 1);
    end

    """Load in the relevant toolboxes for the MLTOOLS.
    
    Description:
    
    """
        importLatest('datasets');
    importLatest('optimi');
    importLatest('ndlutil');
    importLatest('netlab');
    importLatest('mvu');
    importLatest('lmvu');
    importLatest('lle');
    importLatest('jdqr');
    importLatest('isomap');

def model = modelAddDynamics(model, type, varargin):

    """Add a dynamics kernel to the model.
    
    Description:
    
    modelAddDynamics(model, type, ...) adds a dynamics model to a
     model.
     Arguments:
      model - the model to add dynamics to.
      type - the type of dynamics model to add in.
      ... - additional arguments to be passed on creation of the
       dynamics model.
        

    See also
    MODELCREATE


    Copyright (c) 2005, 2006, 2007 Neil D. Lawrence
    
    """
        
    type = [type 'Dynamics'];
    model.dynamics = modelCreate(type, model.q, model.q, model.X, varargin{:});
    params = modelExtractParam(model);
    model = modelExpandParam(model, params);
    

def model = modelCreate(type, numIn, numOut, varargin):

    """Create a model of the specified type.
    
    Description:
    
    model = modelCreate(type, numIn, numOut, P3,...) creates a model
     of the given type.
     Returns:
      model - the model created.
     Arguments:
      type - the type of the model to create, for example, 'kbr' for
       kernel based regression, 'mlp' for multi-layer perceptron,
       'linear' for a linear model.
      numIn - number of inputs to the model (or latent dimensions for
       latent variable models.
      numOut - number of outputs from the model (or data dimensions for
       latent variable models.
      P3,... - optional arguments to be passed to the model creation
       code.
        

    See also
    MODELEXPANDPARAM, MODELEXTRACTPARAM


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([type, 'Create']);
    model = fhandle(numIn, numOut, varargin{:});
def modelDisplay(model, varargin):

    """Display a text output of a model.
    
    Description:
    def modelDisplay(model, varargin):
%
    """
        
    % Check if the model has display code.
    if exist([model.type 'Display'])==2
      fhandle = str2func([model.type 'Display']);
      fhandle(model, varargin{:});
    end

def model = modelExpandParam(model, params, dim):

    """Update a model structure with parameters.
    
    Description:
    
    model = modelExpandParam(model, param) returns a model structure
     filled with the parameters in the given vector. This is used as a
     helper function to enable parameters to be optimised in, for
     example, the NETLAB optimisation functions.
     Returns:
      model - model structure with the given parameters in the relevant
       locations.
     Arguments:
      model - the model structure in which the parameters are to be
       placed.
      param - vector of parameters which are to be placed in the model
       structure.
        
        

    See also
    MODELEXTRACTPARAM, SCG, CONJGRAD


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    
    With modifications by Cark Henrik Ek 2007
    
    """
        
    if isfield(model, 'paramGroups')
      params = params*model.paramGroups';
    end
    
    fhandle = str2func([model.type 'ExpandParam']);
    if(nargin<3)
      model = fhandle(model, params);
    else
      model = fhandle(model,params,dim);
    end
def [params, names] = modelExtractParam(model, dim):

    """Extract the parameters of a model.
    
    Description:
    
    param = modelExtractParam(model) Extract parameters from the model
     into a vector of parameters for optimisation.
     Returns:
      param - vector of parameters extracted from the model.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
        
        
        

    See also
    MODELEXPANDPARAM, SCG, CONJGRAD


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    
    With modifications by Cark Henrik Ek 2007
    
    
    With modifications by Mauricio Alvarez 2008
    
    """
        
    fhandle = str2func([model.type 'ExtractParam']);
    if nargout < 2
        if(nargin<2)
            params = fhandle(model);
        else
            params = fhandle(model,dim);
        end
    else
        [params, namesTemp] = fhandle(model);
    end
    if isfield(model, 'paramGroups')
        paramGroups = model.paramGroups;
        for i = 1:size(paramGroups, 2)
            ind = find(paramGroups(:, i));
            if nargout > 1
                names{i} = namesTemp{ind(1)};
                if length(ind) > 1
                    for j = 2:length(ind)
                        names{i} = [names{i} ', ' namesTemp{ind(j)}];
                    end
                end
            end
            paramGroups(ind(2:end), i) = 0;
        end
        params = params*paramGroups;
    else
        if nargout>1
            names = namesTemp;
        end
    end
    end

def [W, b] = modelGetOutputWeights(model):

    """Wrapper function to return output weight and bias matrices.
    
    Description:
    
    [W, b] = modelGetOutputWeights(model) returns the output weight
     and bias matrices for any mapping model that can be interpreted as
     a generalised linear model (e.g. rbf networks, kernel based
     regressions, multi layer perceptrons, linear).
     Returns:
      W - the output weight matrix.
      b - the output biases.
     Arguments:
      model - the mapping model.
        

    See also
    MLPCREATE, RBFCREATE, KBRCREATE, LINEARCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    switch model.type
     case 'mlp'
      W = model.w2;
      b = model.b2;
     case 'rbf'
      W = model.w2;
      b = model.b2;
     case 'kbr'
      W = model.A;
      b = model.bias;
     case 'linear'
      W = model.W;
      b = model.b;
     otherwise 
      error('Model has no implementation of output weights and biases.');
    end

def g = modelGradient(params, model, varargin):

    """Gradient of error function to minimise for given model.
    
    Description:
    
    g = modelGradient(params, model, ...) gives the gradient of the
     objective function for a model. By default the objective function
     is a negative log likelihood.
     Returns:
      g - the gradient of the error function to be minimised.
     Arguments:
      params - parameter vector to evaluate at.
      model - model structure to optimise.
      ... - optional additional arguments.
        

    See also
    MODELLOGLIKEGRADIENT, MODELOBJECTIVE, MODELOPTIMISE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = [model.type 'Gradient'];
    if exist(fhandle) == 2
      fhandle = str2func(fhandle);
      g = fhandle(params, model, varargin{:});
    else
      fhandle = str2func([model.type 'ExpandParam']);
      model = fhandle(model, params);
      fhandle = str2func([model.type 'LogLikeGradients']);
      g = - fhandle(model, varargin{:});
    end

def modelGradientCheck(model, varargin):

    """Check gradients of given model.
    
    Description:
    
    modelGradientCheck(model, ...) checks the supplied gradient
     function and the supplied objective function to ensure that the
     numerical gradients (as computed with the objective function)
     match the analytically computed gradients.
     Arguments:
      model - the model for which gradients are to be checked.
      ... - additional arguments that are passed to the objective and
       gradient functions (after the parameter vector which is always
       assumed to be the first argument passed).
        

    See also
    MODELOBJECTIVE, MODELGRADIENT, MODELCREATE


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    
    [params, names] = modelExtractParam(model);
    if length(names) == 0
      for i = 1:model.numParams
        names{i} = ['Param ' num2str(i)];
      end
    end
    
    if length(names) ~= length(params)
      error('Names array does not match length of params array');
    end
    
    L = 0;
    change = 1e-6;
    origParams = params;
    for i = 1:length(params)
      params(i) = origParams(i) + change;
      Lplus = modelObjective(params, model, varargin{:});
      params(i) = origParams(i) - change;
      Lminus = modelObjective(params, model, varargin{:});
      diff(i) = (Lplus - Lminus)/(2*change);
      params(i) = origParams(i);
    end
    
    anal = modelGradient(origParams, model, varargin{:});
    
    delta = anal-diff;
    
    
    
    paramMaxDiff = max(max(abs(diff-anal)));
    if paramMaxDiff > 100*change
      l = 0;
      for i = 1:length(names)
        if l < length(names{i})
          l = length(names{i});
        end
      end
      
      fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
      for i = 1:length(names)
        if(abs(delta(i)/max([abs(anal(i)) 1]))>=1e-4)
    
          spaceLen = l - length(names{i});
          space = char(repmat(32, 1, spaceLen));
          fprintf([space names{i} ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
                  anal(i), diff(i), diff(i) - anal(i));
        end
      end
    
    end
    fprintf('Param max diff: %2.6f.\n', paramMaxDiff);
    
    
    

def H = modelHessian(params, model, varargin):

    """Hessian of error function to minimise for given model.
    
    Description:
    
    H = modelHessian(params, model, ...) gives the Hessian of the
     objective function for a model. By default the objective function
     is a negative log likelihood.
     Returns:
      H - the Hessian of the error function to be minimised.
     Arguments:
      params - parameter vector to evaluate at.
      model - model structure to optimise.
      ... - optional additional arguments.
        

    See also
    MODELLOGLIKEHESSIAN, MODELOBJECTIVE, MODELOPTIMISE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = [model.type 'Hessian'];
    if exist(fhandle) == 2
      fhandle = str2func(fhandle);
      H = fhandle(params, model, varargin{:});
    else
      fhandle = str2func([model.type 'ExpandParam']);
      model = fhandle(model, params);
      fhandle = str2func([model.type 'LogLikeHessian']);
      H = - fhandle(model, varargin{:});
    end

def g = modelLogLikeGradients(model):

    """Compute a model's gradients wrt log likelihood.
    
    Description:
    
    g = modelLogLikeGradients(model) is a wrapper function to compute
     the gradients of the log likelihood of a given model.
     Returns:
      g - teh gradients of the likelihood with respect to the
       parameters.
     Arguments:
      model - the model for which likelihoods are computed.
        

    See also
    MODELCREATE


    Copyright (c) 2006, 2005 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.type 'LogLikeGradients']);
    g = fhandle(model);
    
    if isfield(model, 'paramGroups')
      g = g*model.paramGroups;
    end

def ll = modelLogLikelihood(model):

    """Compute a model log likelihood.
    
    Description:
    
    ll = modelLogLikelihood(model) computes the log likelihood of the
     given model.
     Returns:
      ll - the log likelihood of the data given the model.
     Arguments:
      model - the model for which the log likelihood is to be computed.
        

    See also
    MODELLOGLIKEGRADIENTS, MODELCREATE


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.type 'LogLikelihood']);
    ll = fhandle(model);
def err = modelObjective(params, model, varargin):

    """Objective function to minimise for given model.
    
    Description:
    
    err = modelObjective(params, model, ...) gives the objective
     function for a model. By default it is the negative log
     likelihood.
     Returns:
      err - the error function to be minimised.
     Arguments:
      params - parameter vector to evaluate at.
      model - model structure to optimise.
      ... - optional additional arguments.
        

    See also
    MODELLOGLIKELIHOOD, MODELGRADIENT, MODELOPTIMISE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    
    fhandle = [model.type 'Objective'];
    if exist(fhandle) == 2
      fhandle = str2func(fhandle);
      err = fhandle(params, model, varargin{:});
    else
      fhandle = str2func([model.type 'ExpandParam']);
      model = fhandle(model, params);
      fhandle = str2func([model.type 'LogLikelihood']);
      err = - fhandle(model, varargin{:});
    end

def model = modelOptimise(model, varargin):

    """Optimise the given model.
    
    Description:
    
    model = modelOptimise(model, P3...) is a wrapper function that
     optimises a given model.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
      P3... - optional additional arguments.
        
        

    See also
    MODELOBJECTIVE, MODELGRADIENT


    Copyright (c) 2006 Neil D. Lawrence
    
    
    With modifications by Carl Henrik Ek 2007
    
    """
        
    if nargin < 2
      varargin = {}
    end
    fhandle = [model.type 'Optimise'];
    if exist(fhandle)==2
      fhandle = str2func(fhandle);
      model = fhandle(model, varargin{:});
    else
      if ~isfield(model, 'display')
        if length(varargin)< 3
          display = 1;
        else
          display = varargin{3};
        end    
        if length(varargin)<4
          iters = 500;
        else
          iters = varargin{4};
        end
        if length(varargin)<2 | isempty(varargin{2})
        else
          model.y = varargin{2};
        end
        if length(varargin)<1 | isempty(varargin{1})
        else
          model.X = varargin{1};
        end
      end
    
      options = optOptions;
      options(14) = iters;
      options(9) = 1;
      options(1) = display;
      
      
      
      params = modelExtractParam(model);
      if(~isempty(params))
        if isfield(model, 'optimiser')
          optim = str2func(model.optimiser);
        else
          optim = str2func('conjgrad');
        end
        
        params = optim('modelObjective', params,  options, ...
                       'modelGradient', model);
        
        model = modelExpandParam(model, params);
      else
        warning('This Model Has No Parameters To Optimise');
      end
    end
    
    

def options = modelOptions(modelType, varargin):

    """Returns a default options structure for the given model.
    
    Description:
    
    options = modelOptions(modelType, ...) returns a default options
     structure for the given model.
     Returns:
      options - options structure.
     Arguments:
      modelType - the type of model.
      ... - optional additional arguments (dependent on model type).
        

    See also
    MODELCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([modelType 'Options']);
    options = fhandle(varargin{:});
def [Y, Phi] = modelOut(model, X, varargin):

    """Give the output of a model for given X.
    
    Description:
    
    Y = modelOut(model, X) gives the output of the model for a given
     input X. For latent variable models it gives a position in data
     space given a position in latent space.
     Returns:
      Y - output location(s) corresponding to given input locations.
     Arguments:
      model - structure specifying the model.
      X - input location(s) for which output is to be computed.
    
    [Phi, Y] = modelOut(model, X) gives the output of the model for a
     given input X. For latent variable models it gives a position in
     data space given a position in latent space.
     Returns:
      Phi - output basis function(s) corresponding to given input
      Y - output location(s) corresponding to given input locations.
     Arguments:
      model - structure specifying the model.
      X - input location(s) for which output is to be computed.
        
        

    See also
    MODELCREATE


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    
    With modifications by Cark Henrik Ek 2008
    
    """
        
    fhandle = str2func([model.type 'Out']);
    if nargout > 1
      [Y, Phi] = fhandle(model, X, varargin{:});
    else
      Y = fhandle(model, X, varargin{:});
    end
    if(isfield(model,'indexOut')&&~isempty(model.indexOut))
      Y(:,setdiff(1:1:size(Y,2),model.indexOut)) = NaN;
    end

def g = modelOutputGrad(model, X, dim):

    """Compute derivatives with respect to params of model outputs.
    
    Description:
    
    g = modelOutputGrad(model, X) gives the gradients of the outputs
     from the model with respect to the parameters for a given set of
     inputs.
     Returns:
      g - gradients of the model output with respect to the model
       parameters for the given input locations. The size of the returned
       matrix is of dimension number of data x number of parameters x
       number of model outputs (which maintains compatability with
       NETLAB).
     Arguments:
      model - the model structure for which gradients are computed.
      X - input locations where gradients are to be computed.
    
    g = modelOutputGrad(model, X, dim) gives the gradients of the
     outputs from the model with respect to the parameters for a given
     set of inputs.
     Returns:
      g - gradients of the model output with respect to the model
       parameters for the given input locations. The size of the returned
       matrix is of dimension number of data x number of parameters.
     Arguments:
      model - the model structure for which gradients are computed.
      X - input locations where gradients are to be comxfputed.
      dim - the dimension of the model for which gradients are required.
        
        

    See also
    MODELCREATE, MODELLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS, MLPDERIV


    Copyright (c) 2005, 2006 Neil D. Lawrence
    
    
    With modifications by Cark Henrik Ek 2007
    
    """
        
    if(nargin>2)
      fhandle = str2func([model.type 'OutputGrad']);
      g = fhandle(model, X, dim);
    else
      fhandle = str2func([model.type 'OutputGrad']);
      gtemp = fhandle(model, X);
      
      if isfield(model, 'paramGroups')
        g = zeros(size(X, 1), size(model.paramGroups, 2), size(gtemp, 3));
        for i = 1:size(gtemp, 3)
          g = gtemp(:, :, i)*model.paramGroups;
        end
      else 
        g = gtemp;
      end
    end
def g = modelOutputGradX(model, X):

    """Compute derivatives with respect to model inputs of model outputs.
    
    Description:
    
    g = modelOutputGradX(model, X) gives the gradients of the outputs
     from the model with respect to the inputs.
     Returns:
      g - gradients of the model output with respect to the model
       parameters for the given input locations.
     Arguments:
      model - the model structure for which gradients are computed.
      X - input locations where gradients are to be computed.
        

    See also
    MODELCREATE, MODELOUTPUTGRAD, MODELLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.type 'OutputGradX']);
    g = fhandle(model, X);
    

def model = modelParamInit(model, varargin):

    """Initialise the parameters of the model.
    
    Description:
    
    model = modelParamInit(model, ...) initialises the parameters of
     the model with some sensible values.
     Returns:
      model - model with parameters initialised.
     Arguments:
      model - model for which initialisation will be performed.
      ... - optional additional arguments.
        

    See also
    MODELCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.type 'ParamInit']);
    options = fhandle(model, varargin{:});
    

def ll = modelPointLogLikelihood(model, varargin):

    """Compute the log likelihood of a given point.
    
    Description:
    
    ll = modelPointLogLikelihood(model, ...) computes the log
     likelihood of the given model.
     Returns:
      ll - the log likelihood of the given data point.
     Arguments:
      model - the model for which the log likelihood is to be computed.
      ... - additional arguments as required.
        

    See also
    MODELLOGLIKELIHOOD, MODELCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.type 'PointLogLikelihood']);
    ll = fhandle(model, varargin{:});
def [model, lbls] = modelReadFromFID(FID, varargin):

    """Load from a FID produced by C++ code.
    
    Description:
    
    model = modelReadFromFID(FID) loads in from a file stream the data
     format produced by C++ code.
     Returns:
      model - the model loaded in from the file.
     Arguments:
      FID - the file ID from where the data is loaded.
        

    See also
    MODELREADFROMFILE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
      version = readVersionFromFID(FID);
      if version < 0.2
        error('Incorrect file version.')
      end
      
      modelType = readStringFromFID(FID, 'baseType');
      feval = str2func([modelType 'ReadFromFID']);
      model = feval(FID, varargin{:});

def [model, lbls] = modelReadFromFile(fileName, varargin):

    """Read model from a file FID produced by the C++ implementation.
    
    Description:
    
    model = modelReadFromFile(fileName) loads in from a file a model
     produced by C++ code. C++ GP implementation.
     Returns:
      model - the model loaded in from the file.
     Arguments:
      fileName - the file ID from where the data is loaded.
        

    See also
    MODELREADFROMFILE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    
    FID = fopen(fileName);
    if FID==-1
      error(['Cannot find file ' fileName])
    end
    model = modelReadFromFID(FID, varargin{:});
    fclose(FID);

def Y = modelSamp(model, X):

    """Give a sample from a model for given X.
    
    Description:
    def Y = modelSamp(model, X):
%
    """
        
    fhandle = str2func([model.type 'Samp']);
    Y = fhandle(model, X);
def model  = modelSetOutputWeights(model, W, b):

    """Wrapper function to return set output weight and bias matrices.
    
    Description:
    
    model = modelSetOutputWeights(model, W, b) sets the output weight
     and bias matrices for any mapping model that can be interpreted as
     a generalised linear model (e.g. rbf networks, kernel based
     regressions, multi layer perceptrons, linear).
     Returns:
      model - the model with updated weights and bias matrices.
     Arguments:
      model - the mapping model.
      W - the output weight matrix.
      b - the output biases.
        

    See also
    MLPCREATE, RBFCREATE, KBRCREATE, LINEARCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    switch model.type
     case 'mlp'
      model.w2 = W;
      model.b2 = b;
     case 'rbf'
      model.w2 = W;
      model.b2 = b;
     case 'kbr'
      model.A = W;
      model.bias = b;
     case 'linear'
      model.W = W;
      model.b = b;
     otherwise 
      error('Model has no implementation of output weights and biases.');
    end

def modelRet = modelTest(modelType, numIn, varargin):

    """Run some tests on the specified model.
    
    Description:
    
    model = modelTest(modelType) runs some tests on the specified
     model to ensure it is correctly implemented.
     Returns:
      model - the model that was generated for the tests.
     Arguments:
      modelType - type of model to test. For example, 'linear' or 'mlp'.
        

    See also
    MODELCREATE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    %if exist([modelType 'Test']) == 2
    %  feval([modelType 'Test'])
    %end
    if ~isstruct(modelType)
      if nargin < 2
        numIn = 4;
      end
      numData = 20;
      numOut = 3;
      
      % Generate some x positions.
      x = randn(numData, numIn);
      x2 = randn(numData, numIn);
      options = modelOptions(modelType);
      model = modelCreate(modelType, numIn, numOut, options);
      model = modelParamInit(model);
      model.X = x;
      model.y = modelOut(model, x2);
      
      % Set the parameters randomly.
      params = modelExtractParam(model);
      params = randn(size(params))./sqrt(randn(size(params)).^2);
      model = modelExpandParam(model, params);
    else
      model = modelType;
    end
    
    % Check model parameter gradient.
    modelGradientCheck(model, varargin{:});
    
    if exist([model.type 'OutputGrad'])==2;
      epsilon = 1e-6;
      [params, names] = modelExtractParam(model);
      if length(names) == 0
        for i = 1:model.numParams
          names{i} = ['Param ' num2str(i)];
        end
      end
      origParams = params;
      Lplus = zeros(size(model.X, 1), model.numParams, model.outputDim);
      Lminus = zeros(size(model.X, 1), model.numParams, model.outputDim);
      for i = 1:length(params);
        params = origParams;
        params(i) = origParams(i) + epsilon;
        model = modelExpandParam(model, params);
        Lplus(:, i, :) = reshape(modelOut(model, model.X), ...
                                 size(model.X, 1), 1, model.outputDim);
        params(i) = origParams(i) - epsilon;
        model = modelExpandParam(model, params);
        Lminus(:, i, :) = reshape(modelOut(model, model.X), ...
                                  size(model.X, 1), 1, model.outputDim);
      end
      params = origParams;
      model = modelExpandParam(model, params);
      gLDiff = .5*(Lplus - Lminus)/epsilon;
      g = modelOutputGrad(model, model.X);
      
      
      outParamMaxDiff = max(max(max(abs(gLDiff-g))));
      if outParamMaxDiff > 2*epsilon
        l = 0;
        for i = 1:length(names)
          if l < length(names{i})
            l = length(names{i});
          end
        end
        
        fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
        for i = 1:length(names)
          spaceLen = l - length(names{i});
          space = char(repmat(32, 1, spaceLen));
          fprintf([space names{i} ':\t%4.6f\t%4.6f\t%4.6f\n'], ...
                  g(i), gLDiff(i), gLDiff(i) - g(i));
        end
      end
      fprintf('Output param max diff: %2.6f.\n', outParamMaxDiff);
    else
      fprintf('No grad of output with respect to param implemented.\n');
    end
    
    
    if exist([model.type 'OutputGradX'])==2;
      epsilon = 1e-6;
      X = model.X;
      origX = X;
      Lplus = zeros(size(model.X, 1), size(model.X, 2), model.outputDim);
      Lminus = zeros(size(model.X, 1), size(model.X, 2), model.outputDim);
      for i = 1:size(model.X, 2)
        X = origX;
        X(:, i) = origX(:, i) + epsilon;
        Lplus(:, i, :) = reshape(modelOut(model, X), ...
                                 size(model.X, 1), 1, model.outputDim);
        X(:, i) = origX(:, i) - epsilon;
        Lminus(:, i, :) = reshape(modelOut(model, X), ...
                                  size(model.X, 1), 1, model.outputDim);
      end
      X = origX;
      gLDiff = .5*(Lplus - Lminus)/epsilon;
      g = modelOutputGradX(model, X);
      
      
      outputXMaxDiff = max(max(max(abs(gLDiff-g))));
      if outputXMaxDiff > 2*epsilon
        fprintf('gX\n')
        disp(g)
        fprintf('gXDiff\n')
        disp(gLDiff)
      end
      fprintf('X max diff: %2.6f.\n', outputXMaxDiff);
    else
      fprintf('No grad of output with respect to X implemented.\n');
    end
    if nargout > 0
      modelRet = model;
    else
      modelDisplay(model);
    end
    

def model = modelTieParam(model, paramsList):

    """Tie parameters of a model together.
    
    Description:
    
    model = modelTieParam(model, paramsList) groups of parameters of a
     model to be seen as one parameter during optimisation of the
     model.
     Returns:
      model - the model with the parameters grouped together.
     Arguments:
      model - the model for which parameters are being tied together.
      paramsList - indices of parameteres to group together. The indices
       are provided in a cell array. Each cell in the array contains a
       vector of indices of parameters that should be considered as one
       parameter. Each group of parameters in each cell should obviously
       be mutually exclusive.
        

    See also
    MODELEXTRACTPARAM, MODEEXPANDPARAM, MODELLOGLIKEGRADIENTS


    Copyright (c) 2003, 2006, 2008 Neil D. Lawrence
    
    """
        
      if ~isfield(model, 'paramGroups')
        if isfield(model, 'nParams')
          model.paramGroups = speye(model.nParams);
        elseif isfield(model, 'numParams')
          model.paramGroups = speye(model.numParams);
        else
          error('Model does not list number of parameters.');
        end
      end
      colToDelete = [];
      for i = 1:length(paramsList)
        
        paramIndices=sort(paramsList{i});
        if any(paramIndices(1) == colToDelete)
          error('Parameter is already being tied')
        end
        for j = 2:length(paramIndices)
          
          model.paramGroups(paramIndices(j), paramIndices(1)) = 1;
          if any(paramIndices(j) == colToDelete)
            error('Parameter has already been tied')
          end
          colToDelete = [colToDelete paramIndices(j)];
        end
      end
      
      model.paramGroups(:, colToDelete) = [];
      if isfield(model, 'nParams')
        % Update to the new number of parameters.
        model.nParams = size(model.paramGroups, 2);
      elseif isfield(model, 'numParams')
        model.numParams = size(model.paramGroups, 2);
      end
    end
def modelWriteToFID(FID, model):

    """Write to a stream a given model.
    
    Description:
    
    modelWriteToFID(FID, model) loads in from a file stream the data
     format produced by C++ code.
     Arguments:
      FID - the file ID from where the data is loaded.
      model - the model loaded in from the file.
        

    See also
    MODELREADFROMFID


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    writeVersionToFID(FID, 0.2);
    modelType = readStringFromFID(FID, 'baseType');
    feval = str2func([modelType 'ReadFromFID']);
    model = feval(FID, varargin{:});

def model = mogCreate(latentDim, dataDim, Y, options):

    """Create a mixtures of Gaussians model.
    
    Description:
    
    model = mogCreate(latentDim, dataDim, Y, options) creates a
     mixtures of probabilistic PCA model.
     Returns:
      model - the initialised mixtures of probabilistic PCA model.
     Arguments:
      latentDim - the latent dimensionality of the components of the
       probabilistic PCA model.
      dataDim - the dimensionality of the data.
      Y - the data to be modelled.
      options - options structure containing the default options.
        

    See also
    GMM, MODELCREATE


    Copyright (c) 2006, 2008 Neil D. Lawrence
    
    """
        
    model.type = 'mog';
    model.covtype = options.covtype;
    model.q = latentDim;
    model.d = dataDim;
    model.N = size(Y, 1);
    
    model.m = options.numComponents;
    model.isInfinite = options.isInfinite;
    if model.isInfinite
      model.a0 = 1;
      model.a1 = options.a1;
    end
    ind = randperm(model.N);
    ind = ind(1:model.m);
    model.Y = Y;
    
    % Initialise means at random and posteriors heuristically.
    model.mean = model.Y(ind, :);
    model.prior = repmat(1/model.m, 1, model.m);
    
    dists = dist2(model.Y, model.mean);
    model.posterior = dists./repmat(sum(dists, 2), 1, model.m);
    
    model = mogUpdateCovariance(model);
def model = mogEstep(model):

    """Do an E-step on an MOG model.
    
    Description:
    
    model = mogEstep(model) carries out an expectation step on a
     mixtures of Gaussians model.
     Returns:
      model - the model with updated posteriors.
     Arguments:
      model - the model which is to be updated.
        

    See also
    MOGCREATE, MOGUPDATEMEAN, MOGUPDATECOVARIANCE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    for i = 1:model.m
      centredY = model.Y - repmat(model.mean(i, :), model.N, 1);
      switch model.covtype
       case 'ppca'
        centredY = centredY/model.U{i};
        logDetTerm = -0.5*logdet([], model.U{i});
        
       case 'spherical'
        centredY = centredY/sqrt(model.sigma2(i));
        logDetTerm = -0.5*model.d*log(model.sigma2(i));
      end
      centredY = sum(centredY.*centredY, 2);
      model.posterior(:, i) = log(model.prior(i))+logDetTerm-0.5* ...
          centredY;
    end
    maxL = max(model.posterior, [], 2);
    model.posterior = model.posterior - repmat(maxL, 1, model.m);
    model.lnposterior = model.posterior;
    model.posterior = exp(model.posterior);
    sumPosterior = repmat(sum(model.posterior, 2), 1, model.m);
    model.posterior = model.posterior./sumPosterior;
    model.lnposterior = model.lnposterior - log(sumPosterior);

def ll = mogLogLikelihood(model):

    """Mixture of Gaussian's log likelihood.
    
    Description:
    
    lll = mogLogLikelihood(model) computes the variational lower bound
     on the log likelihood of a mixtures of probabilistic PCA model, it
     wraps the mogLowerBound command.
     Returns:
      lll - the lower bound on the log likelihood computed for the
       model.
     Arguments:
      model - the model for which log likelihood is to be computed.
        

    See also
    MOGCREATE, MOGLOWERBOUND


    Copyright (c) 2006, 2008 Neil D. Lawrence
    
    """
        model = mogEstep(model);
    ll = mogLowerBound(model);
    

def lll = mogLowerBound(model):

    """Computes lower bound on log likelihood for an MOG model.
    
    Description:
    
    lll = mogLowerBound(model) computes the variational lower bound on
     the log likelihood of a mixtures of probabilistic PCA model.
     Returns:
      lll - the lower bound on the log likelihood computed for the
       model.
     Arguments:
      model - the model for which log likelihood is to be computed.
        

    See also
    MOGCREATE, MODELLOGLIKELIHOOD


    Copyright (c) 2006, 2008 Neil D. Lawrence
    
    """
        
    lll = -sum(sum(xlogy(model.posterior)));
    
    if model.isInfinite
      % DP add in extra posterior entropy etc.
    end
    
    for i = 1:model.m
      centredY = model.Y - repmat(model.mean(i, :), model.N, 1);
      switch model.covtype
       case 'ppca'
        centredY = centredY/model.U{i};
        logDetTerm = logdet([], model.U{i});
       case 'spherical'
        centredY = centredY/sqrt(model.sigma2(i));
        logDetTerm = model.d*log(model.sigma2(i));
      end
      centredY = sum(centredY.*centredY, 2);
      lll = lll - 0.5*sum(model.posterior(:, i).*(logDetTerm+centredY-2*log(model.prior(i)+1e-300)));
    end
    lll = lll - model.N*model.d/2*log(2*pi);  

def [m, C] = mogMeanCov(model):

    """Project a mixture of Gaussians to a low dimensional space.
    
    Description:
    
    [m, C] = mogMeanCov(model) computes the mean and covariance of the
     distribution given by a mixture of Gaussians.
     Returns:
      m - mean of the mixture of Gaussians.
      C - covariance of the mixture of Gaussians.
     Arguments:
      model - model for which mean and covariance is required.
        

    See also
    MOGCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    m = zeros(1, model.d);
    modelSecondMoment = zeros(model.d);
    for i = 1:model.m
      m = m + model.prior(i)*model.mean(i, :);
    end
    switch model.covtype
      case 'ppca'
       for i = 1:model.m
         modelSecondMoment = modelSecondMoment ...
             + model.prior(i)*model.mean(i, :)'*model.mean(i, :) ...
             + model.prior(i)*(model.W{i}*model.W{i}' + eye(model.d)* ...
                               model.sigma2(i));
       end
      case 'spherical'
       for i = 1:model.m
         modelSecondMoment = modelSecondMoment ...
             + model.prior(i)*model.mean(i, :)'*model.mean(i, :) ...
             + model.prior(i)*(eye(model.d)*model.sigma2(i));
       end
    end
    C = modelSecondMoment - (m'*m);
def model = mogOptimise(model, display, iters):

    """Optimise an MOG model.
    
    Description:
    
    model = mogOptimise(model) optimises an mixtures of Gaussians
     model via the expectation maximisation algorithm.
     Returns:
      model - the optimised model.
     Arguments:
      model - the model to be optimised.
        

    See also
    MMPCACREATE, MODELOPTIMISE


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    diffll = 1;
    iter = 0;
    ll = mogLowerBound(model);
    while abs(diffll)>1e-6 & iter<iters
      iter = iter + 1;
      model =  mogEstep(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'E-step');
      end
      model = mogUpdatePrior(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'Prior Update');
      end
      model = mogUpdateMean(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'Mean Update');
      end
      model = mogUpdateCovariance(model);
      if display > 1
        [ll, oldll] = boundCheck(model, ll, 'Covariance Update');
      end
      if display > 1
      else
        oldll = ll;
        ll = mogLowerBound(model);
      end
      diffll = ll -oldll;
      if display
        fprintf('Iteration %d log-likelihood: %2.6f\n', iter, ll)
      end
    end
    
    
def [ll, oldll] = boundCheck(model, oldll, step):
    
    % BOUNDCHECK Helper function for checking bound.
    
    ll = mogLowerBound(model);
    diffll = ll - oldll;
    if ll -oldll < 0
      warning(['Log likelihood went down by ' num2str(diffll) 'in ' ...
               'step: ' step ' in mogOptimise'])
    end

def options = mogOptions(numComp):

    """Sets the default options structure for MOG models.
    
    Description:
    
    options = mogOptions(numComponents) sets the default options
     structure for mixtures of Gaussians models.
     Returns:
      options - structure containing the default options.
     Arguments:
      numComponents - number of components in the mixture model.
        

    See also
    MOGCREATE


    Copyright (c) 2006, 2008 Neil D. Lawrence
    
    """
        
    options.numComponents = numComp;
    options.covtype = 'ppca';
    % Whether it is an infinite mixture (false by default);
    options.isInfinite = false;
    options.a1=1;
def mogPrintPlot(model, lbls, capName, experimentNo):

    """Print projection of MOG into two dimensions.
    
    Description:
    
    mogPrintPlot(model, lbls, capName, experimentNo) prints a
     projection of mixtures of Gaussians into two dimensions.
     Arguments:
      model - the model to use for plotting the latent space.
      lbls - any lables that are available for plotting.
      capName - the name of the saved plots.
      experimentNo - the experiment number to assign to the files.
        

    See also
    MOGSCATTERPLOT


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    if model.d>2
      model = mogProject(model, 2);
    end
    
    modelType = model.type;
    modelType(1) = upper(modelType(1));
    
    
    fileName = ['dem' capName modelType num2str(experimentNo)];
    
    clf
    ax = axes('position', [0.05 0.05 0.9 0.9]);
    hold on
    if ~isempty(lbls) && ~strcmp(lbls, 'connect')
      mogTwoDPlot(model, lbls, getSymbols(size(lbls, 2)));
    else
      mogTwoDPlot(model, lbls);
    end
    
    
    piVals = linspace(-pi, pi, 200)';
    for i=1:model.m
      a = line(model.mean(i, 1), model.mean(i, 2), 'marker', 'o');
      set(a, 'linewidth', 2, 'markersize', 10)
      x = [sin(piVals) cos(piVals)];
      el = x*model.U{i};
      line(model.mean(i, 1) + el(:, 1), model.mean(i, 2) + el(:, 2), ...
          'linewidth', 2);
    end
    xLim = [min(model.Y(:, 1)) max(model.Y(:, 1))]*1.1;
    yLim = [min(model.Y(:, 2)) max(model.Y(:, 2))]*1.1;
    set(ax, 'xLim', xLim);
    set(ax, 'yLim', yLim);
    set(gca, 'fontsize', 20);
    print('-depsc', ['../tex/diagrams/' fileName])
    print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
    
    % make smaller for PNG plot.
    pos = get(gcf, 'paperposition')
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    fontsize = get(gca, 'fontsize');
    set(gca, 'fontsize', fontsize/2);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    print('-dpng', ['../html/' fileName])
    set(gcf, 'paperposition', origpos);
    
    figure
    clf
    ax = axes('position', [0.05 0.05 0.9 0.9]);
    hold on
    if ~isempty(lbls) && ~strcmp(lbls, 'connect')
      mogTwoDPlot(model, lbls, getSymbols(size(lbls, 2)));
    else
      mogTwoDPlot(model, lbls);
    end
    
    %xLim = [min(model.Y(:, 1)) max(model.Y(:, 1))]*1.1;
    %yLim = [min(model.Y(:, 2)) max(model.Y(:, 2))]*1.1;
    set(ax, 'xLim', xLim);
    set(ax, 'yLim', yLim);
    
    set(ax, 'fontname', 'arial');
    set(ax, 'fontsize', 20);
    print('-depsc', ['../tex/diagrams/' fileName 'NoOvals'])
    print('-deps', ['../tex/diagrams/' fileName 'NoOvalsNoColour'])

def model = mogProject(model, dimension):

    """Project a mixture of Gaussians to a low dimensional space.
    
    Description:
    
    model = mogProject(model, dimension) projects a mixture of
     Gaussians down to a lower dimensional space (typically for
     visualisation).
     Returns:
      model - the reduced dimensional mixture of Gaussians.
     Arguments:
      model - the mixture of Gaussians to project down.
      dimension - the dimension to project to.
        

    See also
    MOGCREATE, MOGTWODPLOT


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
      
    [m, C] = mogMeanCov(model);
    [U, V] = eig(C);
    
    v = diag(V);
    [v, ind] = sort(v);
    ind = ind(end:-1:1);
    U = U(:, ind);
    
    Uq = U(:, 1:dimension);
    
    model.Y = model.Y*Uq;
    model.mean = model.mean*Uq;
    model.d = dimension;
    
    switch model.covtype
     case 'ppca'
      for i = 1:model.m
        model.W{i} = (model.W{i}'*Uq)';
        model.U{i} = sqrt(model.sigma2(i))*eye(model.d);
        for j = 1:model.q
          model.U{i} = cholupdate(model.U{i}, model.W{i}(:, j));
        end
      end
    
    end

def x = mogSample(model, numSamples):

    """Sample from a mixture of Gaussians model.
    
    Description:
    
    x = mogSample(model, numSamples) samples from a mixture of
     Gaussians.
     Returns:
      x - the samples from the model.
     Arguments:
      model - the model that you want to sample from.
      numSamples - the number of samples required.
        

    See also
    MOGCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
    p = rand(numSamples, 1);
    bins = cumsum(model.prior);
    compNo = sum(repmat(p, 1, model.m)<repmat(bins, numSamples, 1), 2);
    x = zeros(numSamples, model.d);
    
    for i = 1:model.m
      ind = find(compNo == i);
      x(ind, :) = repmat(model.mean(i, :), length(ind), 1);
      switch model.covtype
       case 'ppca'
        samps = randn(length(ind), model.q);
        samps = samps*model.W{i}';
        samps = samps + randn(length(ind), model.d)*sqrt(model.sigma2(i));
        x(ind, :) = x(ind, :) + samps;
       case 'spherical'
        x(ind, :) = x(ind, :) + randn(length(ind), model.d)* ...
            sqrt(model.sigma(i));
      end
    end

def returnVal = mogTwoDPlot(model, lbl, symbol):

    """Helper function for plotting the labels in 2-D.
    
    Description:
    
    mogTwoDPlot(model, lbl, symbol) helper function for plotting an
     MOG in 2-D with symbols.
     Arguments:
      model - the data to plot.
      lbl - the labels of the data point.
      symbol - the symbols to use for the different labels.
        

    See also
    MOGSCATTERPLOT


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    if nargin < 2
      lbl = [];
    end
    if(strcmp(lbl, 'connect'))
      connect = true;
      lbl = [];
    else
      connect = false;
    end
    if model.d>2
      mod2 = mogProject(model, 2);
    else
      mod2 = model;
    end
    
    if nargin < 3
      if isempty(lbl)
        symbol = getSymbols(1);
      else
        symbol = getSymbols(size(lbl,2));
      end
    end
    
    axisHand = gca;
    returnVal = [];
    nextPlot = get(axisHand, 'nextplot');
    if ~isempty(lbl)
      for i = 1:size(mod2.Y, 1)
        if i == 2
          set(axisHand, 'nextplot', 'add');
        end
        labelNo = find(lbl(i, :));
        try 
          returnVal = [returnVal; plot(mod2.Y(i, 1), mod2.Y(i, 2), symbol{labelNo})];
        catch
          if strcmp(lasterr, 'Index exceeds matrix dimensions.')
    	error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
          end
        end
      end
      set(axisHand, 'nextplot', nextPlot);
    else
      if connect
        returnVal = plot(mod2.Y(:, 1), mod2.Y(:, 2), 'rx-');
      else
        returnVal = plot(mod2.Y(:, 1), mod2.Y(:, 2), 'rx');
      end
    end
    set(returnVal, 'markersize', 10);
    set(returnVal, 'linewidth', 2);
    

def model = mogUpdateCovariance(model):

    """Update the covariances of an MOG model.
    
    Description:
    
    model = mogUpdateCovariance(model) updates the covariance matrices
     of a mixtures of Gaussians model. The implementation currently
     uses an eigenvalue based update.
     Returns:
      model - the model with updated covariances.
     Arguments:
      model - the model which is to be updated.
        

    See also
    MOGCREATE, MOGUPDATEMEAN, MOGESTEP


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    for i = 1:model.m
    
      centredY = model.Y - repmat(model.mean(i,:), model.N, 1);
      centredY = centredY.*repmat(sqrt(model.posterior(:,i)), 1, model.d);
      switch model.covtype
        case 'ppca'
         C = (centredY'*centredY+0.001*eye(model.d))/sum(model.posterior(:, i)+.001);
         [vec, val] = eig(C);
         val = diag(val);
         [val, ind] = sort(val);
         ind = ind(end:-1:1);
         val = val(end:-1:1);
         vec = vec(:, ind(1:model.q));
         sigma2 = mean(val(model.q+1:end));
         if sigma2<eps
           sigma2 = eps;
         end
         lambda = val(1:model.q) - sigma2;
         %[sigma2, eigVec, lambda] = ppca(C, model.q);
         if length(lambda) ~= model.q
           % Something's wrong here ...
           sigma2 = 1e-6;
           warning('Not enough eigenvalues extracted.')
           lambdaTemp = lambda;
           lambda = zeros(model.q, 1);
           lambda(1:length(lambdaTemp)) = lambdaTemp;
         end 
        
         model.sigma2(i) = sigma2;
         model.W{i} = vec*diag(sqrt(lambda));
         model.U{i} = sqrt(sigma2)*eye(model.d);
         for j = 1:model.q
           model.U{i} = cholupdate(model.U{i}, model.W{i}(:, j));
         end
       case 'spherical'
        model.sigma2(i) = sum(sum(centredY.*centredY))/(model.d*sum(model.posterior(:, i)));
      end
    end        
    

def model = mogUpdateMean(model):

    """Update the means of an MOG model.
    
    Description:
    
    model = mogUpdateMean(model) updates the mean vectors of a
     mixtures of Gaussians model.
     Returns:
      model - the model with updated means.
     Arguments:
      model - the model which is to be updated.
        

    See also
    MOGCREATE, MOGUPDATECOVARIANCE, MOGESTEP


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    for i = 1:model.m
      if sum(model.posterior(:, i)) ~= 0
        for j = 1:model.d
          model.mean(i, j) = model.posterior(:, i)'*model.Y(:, ...
                                                            j)/sum(model.posterior(:, i));
        end    
      else
        p = exp(model.lnposterior(:, i) - max(model.lnposterior(:, i)));
        for j = 1:model.d
          model.mean(i, j) = p'*model.Y(:, j)/sum(p);
        end
      end
    end
    

def model = mogUpdatePrior(model):

    """Update the priors of an MOG model.
    
    Description:
    
    model = mogUpdatePrior(model) updates the prior probabilities of a
     mixtures of Gaussians model.
     Returns:
      model - the model with updated priors.
     Arguments:
      model - the model which is to be updated.
        

    See also
    MOGCREATE, MOGUPDATEMEAN, MOGUPDATECOVARIANCE, MOGESTEP


    Copyright (c) 2006, 2008 Neil D. Lawrence
    
    """
        
    if model.isInfinite
      % First compute expectations of v.
      sumS = sum(model.posterior);
      a0bar = model.a0 + sumS; % Posterior value for a_0.
      a1bar = model.a1 + cumsum(sumS); % Posterior value for a_1.
      model.v = a0bar./(a0bar+a1bar);
      tmp = cumprod(1-model.v);
      model.prior = model.v;
      model.prior(2:end) = model.prior(2:end).*tmp(1:end-1);
    else
      model.prior = mean(model.posterior);
      model.prior(find(model.prior==0))=1e-100;
    end

def model = multimodelCreate(inputDim, outputDim, varargin):

    """Create a MULTIMODEL model.
    
    Description:
        The MULTIMODEL is a way of performing multi-task learning by sharing
        model parameters across a range of models. The default (simple)
        assumption is that the data is conditionally independent given the
        parameters, i.e. the log likelihood is the sum of the log likelihood of
        the models.
        
        
    
    model = multimodelCreate(inputDim, outputDim, options) creates a
     multi-task learning wrapper model structure given an options
     structure.
     Returns:
      model - the model structure with the default parameters placed in.
     Arguments:
      inputDim - the input dimension of the model.
      outputDim - the output dimension of the model.
      options - an options structure that determines the form of the
       model.
        

    See also
    MODELCREATE, MULTIMODELOPTIONS, MULTIMODELPARAMINIT, MODELCREATE


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        options = varargin{end};
    model.numModels = options.numModels;
    model.type = 'multimodel';
    model.compType = options.type;
    model.inputDim = inputDim;
    model.outputDim = outputDim;
    % Indices of parameters to be trained separately for each model.
    model.separateIndices = options.separate;
    model.numSep = length(model.separateIndices);
    for i = 1:model.numModels
      varargput = cell(1, length(varargin)-1);
      for j = 1:length(varargput)
        varargput{j} = varargin{j}{i};
      end
      model.comp{i} = modelCreate(model.compType, inputDim, outputDim, ...
                                  varargput{:}, options.compOptions);
    end
    if isfield(model.comp{i}, 'numParams');
      model.numParams = model.comp{1}.numParams;
    else
      model.numParams = length(modelExtractParam(model.comp{1}));
    end
    model.sharedIndices = 1:model.numParams;
    model.sharedIndices(model.separateIndices) = [];
    
    model.numParams = model.numParams + (model.numModels-1)*model.numSep;
    
    if isfield(options, 'optimiser') && ~isempty(options.optimiser)
        model.optimiser = options.optimiser;    
    end
    

def multimodelDisplay(model, spacing):

    """Display parameters of the MULTIMODEL model.
    
    Description:
    
    multimodelDisplay(model) displays the parameters of the multi-task
     learning wrapper model and the model type to the console.
     Arguments:
      model - the model to display.
    
    multimodelDisplay(model, spacing)
     Arguments:
      model - the model to display.
      spacing - how many spaces to indent the display of the model by.
        

    See also
    MULTIMODELCREATE, MODELDISPLAY


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Multi-model:\n')
    for i = 1:length(model.comp)
      fprintf('Component %d:\n', i)
      spacing = length(spacing)+2;
      modelDisplay(model.comp{i}, spacing);
    end

def model = multimodelExpandParam(model, params):

    """Create model structure from MULTIMODEL model's parameters.
    
    Description:
    
    model = multimodelExpandParam(model, param) returns a multi-task
     learning wrapper model structure filled with the parameters in the
     given vector. This is used as a helper function to enable
     parameters to be optimised in, for example, the NETLAB
     optimisation functions.
     Returns:
      model - model structure with the given parameters in the relevant
       locations.
     Arguments:
      model - the model structure in which the parameters are to be
       placed.
      param - vector of parameters which are to be placed in the model
       structure.
        
        

    See also
    MULTIMODELCREATE, MULTIMODELEXTRACTPARAM, MODELEXPANDPARAM


    Copyright (c) 2007, 2008 Neil D. Lawrence
    
    
    With modifications by Mauricio Alvarez 2008
    
    """
        
    endVal = model.numParams - model.numModels*model.numSep;
    baseParams = params(1:endVal);
    passParams(model.sharedIndices) = baseParams;
    if ~isempty(model.separateIndices)
        for i = 1:length(model.comp)
            startVal = endVal + 1;
            endVal = endVal + model.numSep;
            passParams(model.separateIndices) = params(startVal:endVal);
            model.comp{i} = modelExpandParam(model.comp{i}, passParams);
        end
    else
        for i = 1:length(model.comp)
            model.comp{i} = modelExpandParam(model.comp{i}, passParams);
        end
    end
    end
def [passParams, passNames] = multimodelExtractParam(model):

    """Extract parameters from the MULTIMODEL model structure.
    
    Description:
    
    param = multimodelExtractParam(model) extracts parameters from the
     multi-task learning wrapper model structure into a vector of
     parameters for optimisation.
     Returns:
      param - vector of parameters extracted from the model.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
    
    [param, names] = multimodelExtractParam(model) extracts parameters
     and parameter names from the multi-task learning wrapper model
     structure.
     Returns:
      param - vector of parameters extracted from the model.
      names - cell array of strings containing names for each parameter.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
        
        
        

    See also
    % SEEALSO MULTIMODELCREATE, MULTIMODELEXPANDPARAM, MODELEXTRACTPARAM, SCG, CONJGRAD


    Copyright (c) 2007, 2008 Neil D. Lawrence
    
    
    With modifications by Mauricio Alvarez 2008, 2009
    
    """
          passParams = zeros(1, model.numParams);
      passNames = cell(1, model.numParams);
      if nargout > 1
        [receiveParams, receiveNames] = modelExtractParam(model.comp{1});
      else
        receiveParams = modelExtractParam(model.comp{1});
      end
      %endVal = model.numParams - model.numSep;
      endVal = model.numParams - model.numSep*model.numModels;
      passParams(1:endVal) = receiveParams(model.sharedIndices);
      if nargout > 1
        % MAURICIO : I think this is wrong. But I didn't change because I'm not
        % sure
        % passNames{1:endVal} = receiveNames{model.sharedIndices};
        passNames(1:endVal) = receiveNames(model.sharedIndices);
      end
      if ~isempty(model.separateIndices)
        startVal = endVal + 1;
        endVal = endVal + model.numSep;
        passParams(startVal:endVal) = receiveParams(model.separateIndices);
        if nargout > 1
          passNames(startVal:endVal) = receiveNames(model.separateIndices);
          for j=startVal:endVal
              passNames{j} = [model.type ' 1 '  passNames{j}];
          end
        end
        for i = 2:length(model.comp)
          startVal = endVal+1;
          endVal = endVal + model.numSep;
          if nargout > 1
            [receiveParams, receiveNames] = modelExtractParam(model.comp{i});
          else
            receiveParams = modelExtractParam(model.comp{i});
          end
          passParams(startVal:endVal) = receiveParams(model.separateIndices);
          if nargout > 1
            passNames(startVal:endVal) = receiveNames(model.separateIndices);
            for j=startVal:endVal
                passNames{j} = [model.type ' ' num2str(i) ' '  passNames{j}];
            end
          end
        end
      end
    end
def g = multimodelLogLikeGradients(model):

    """Gradient of MULTIMODEL model log likelihood with respect to parameters.
    
    Description:
    
    g = multimodelLogLikeGradients(model) computes the gradient of the
     multi-task learning wrapper model's log likelihood with respect to
     the parameters.
     Returns:
      g - the returned gradients.
     Arguments:
      model - model structure for which gradients are being computed.
        
        

    See also
    % SEEALSO MULTIMODELCREATE, MULTIMODELLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS


    Copyright (c) 2007, 2008 Neil D. Lawrence
    
    
    With modifications by Mauricio Alvarez 2009
    
    """
        
      g = zeros(1, model.numParams);
      endShared = model.numParams - model.numModels*model.numSep;
      endVal = endShared; 
      for i = 1:length(model.comp)
        gModel = modelLogLikeGradients(model.comp{i}); 
        g(1:endShared) = g(1:endShared) + gModel(model.sharedIndices);
        if ~isempty(model.separateIndices)
          startVal = endVal + 1;
          endVal = endVal + model.numSep;
          %g(startVal:endVal) = g(model.separateIndices);
          g(startVal:endVal) = gModel(model.separateIndices);      
        end
      end
    end

def ll = multimodelLogLikelihood(model):

    """Log likelihood of MULTIMODEL model.
    
    Description:
    
    ll = multimodelLogLikelihood(model) computes the log likelihood of
      the multi-task learning wrapper model.
     Returns:
      ll - the computed log likelihood.
     Arguments:
      model - the model structure for which log likelihood is being
       computed.
        

    See also
    MULTIMODELCREATE, MULTIMODELLOGLIKEGRADIENTS, MODELLOGLIKELIHOOD


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    ll = 0;
    for i = 1:length(model.comp)
      ll = ll + modelLogLikelihood(model.comp{i});
    end
def options = multimodelOptions(modelType, number, varargin):

    """Create a default options structure for the MULTIMODEL model.
    
    Description:
    
    options = multimodelOptions(modelType, number, ...) creates a
     default options structure for the multi-task learning wrapper
     model structure.
     Returns:
      options - the default options structure.
     Arguments:
      modelType - the model type that the multi-task model is based on.
      number - the number of components in the multi-task model.
      ... - optional additional arguments to be passed to the
       'sub-model's options structure.
        

    See also
    MULTIMODELCREATE, MODELOPTIONS


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    options.type = modelType;
    options.numModels = number;
    fhandle = str2func([modelType 'Options']);
    options.compOptions = fhandle(varargin{:});
def model = multimodelParamInit(model):

    """MULTIMODEL model parameter initialisation.
    
    Description:
    
    model = multimodelParamInit(model) initialises the multi-task
     learning wrapper model structure with some default parameters.
     Returns:
      model - the model structure with the default parameters placed in.
     Arguments:
      model - the model structure which requires initialisation.
        

    See also
    MULTIMODELCREATE, MODELCREATE, MODELPARAMINIT


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    for i = 1:length(model.comp)
      model.comp{i} = modelParamInit(model.comp{i});
    end
def X = mvuEmbed(Y,dims,k):

    """Embed data set with MVU.
    
    Description:
    
    X = mvuEmbed(Y, dims, k) Embed a given data set with Weinberger et
     al.'s MVU algorithm.
     Returns:
      X - embedding
     Arguments:
      Y - Data
      dims - Dimensionality of Embedding (default = 2)
      k - Number of Neighbours in Proximity Graph (default = 7)
        

    See also
    PPCAEMBED, LLEEMBED, LMVUEMBED


    Copyright (c) Neil D. Lawrence, 2007 Carl Henrik Ek
    
    """
        
    if(nargin<3)
      k = 7;
      if(nargin<2)
        dims = 2;
        if(nargin<1)
          error('To Few Arguments');
        end
      end
    end
    
    if(any(any(isnan(Y))))
      error('Cannot run MVU when missing data is present.');
    end
    
    X = mvu(distance(Y'),k);
    
    X = X(1:1:dims,:)';
    
    
def D = distance(Y):
      
      D = sqrt(dist2(Y', Y'));
    return
def ind = paramNameRegularExpressionLookup(model, pattern):

    """Returns the indices of the parameter containing the given regular expression.
    
    Description:
    
    ind = paramNameRegularExpressionLookup(model, pattern) returns the
     index of the parameters which contain the given regular
     expression. If no matches are found then an empty vector is
     returned.
     Returns:
      ind - the indices of those parameters in the model.
     Arguments:
      model - the model for which parameters are reverse looked up.
      pattern - the regular expression that should match the names.
        

    See also
    CMPNDTIEPARAMETERS, PARAMNAMEREVERSELOOKUP


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
      ind = [];
      [void, names] = modelExtractParam(model);
      for i = 1:length(names)
        if(regexp(names{i}, pattern))
          ind = [ind i];
        end
      end
    end
def ind = paramNameReverseLookup(model, name):

    """Returns the index of the parameter with the given name.
    
    Description:
    
    ind = paramNameReverseLookup(model, name) returns the index of the
     parameter with the given name. If no matches are found then an
     empty vector is returned.
     Returns:
      ind - the indices of those parameters in the model.
     Arguments:
      model - the model for which parameters are reverse looked up.
      name - the name of the parameter that is sought.
        

    See also
    CMPNDTIEPARAMETERS, PARAMNAMEREGULAREXPRESSIONLOOKUP


    Copyright (c) 2008 Neil D. Lawrence
    
    """
          
      ind = [];
      [void, names] = modelExtractParam(model);
      for i = 1:length(names)
        if(strcmp(names{i}, name))
          ind = i;
          return
        end
      end
    end
def model = ppcaCreate(inputDim, outputDim, Y, options):

    """Density network model.
    
    Description:
    
    model = ppcaCreate(inputDimension, outputDim, Y, options) creates
     a structure for a density network.
     Returns:
      model - model structure containing the neural network specified.
     Arguments:
      inputDimension - dimension of input data.
      outputDim - dimension of target data.
      Y - the data to be modelled in design matrix format (as many rows
       as there are data points).
      options - options structure as returned by ppcaCreate.
        

    See also
    PPCAOPTIONS, MLPCREATE, RBFCREATE, KBRCREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    
    model.type = 'ppca';
    
    if size(Y, 2) ~= outputDim
      error(['Input matrix Y does not have dimension ' num2str(d)]);
    end
    model.X = ppcaEmbed(Y, inputDim);
    model.y = Y;
    model.q = inputDim;
    model.d = outputDim;
    model.N = size(Y, 1);
    model.b = mean(model.y);
    model.W = pdinv(model.X'*model.X)*model.X'*(model.y-repmat(model.b, model.N, ...
                                                      1));
    yDiff = model.y - (model.X*model.W + repmat(model.b, model.N, 1));
    yDiff = yDiff.*yDiff;
    model.beta = model.N*model.d/sum(sum(yDiff));
def [X, sigma2] = ppcaEmbed(Y, dims):

    """Embed data set with probabilistic PCA.
    
    Description:
    
    [X, sigma2] = ppcaEmbed(Y, dims) returns latent positions for a
     given data set via probabilistic PCA.
     Returns:
      X - the latent positions.
      sigma2 - the variance not explained by the latent positions.
     Arguments:
      Y - the data set which you want the latent positions for.
      dims - the dimensionality of the latent space.
        

    See also
    LLEEMBED, ISOMAPEMBED


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    if ~any(any(isnan(Y)))
      [v, u] = pca(Y);
      v(find(v<0))=0;
      Ymean = mean(Y);
      Ycentre = zeros(size(Y));
      for i = 1:size(Y, 2);
        Ycentre(:, i) = Y(:, i) - Ymean(i);
      end
      X = Ycentre*u(:, 1:dims)*diag(1./sqrt(v(1:dims)));
      sigma2 = mean(v(dims+1:end));
    else
      % Hacky implementation of Probabilistic PCA for when there is missing data.
      iters = 100;
      % Initialise W randomly
      d = size(Y, 2);
      q = dims;
      N = size(Y, 1);
      W = randn(d, q)*1e-3;
      sigma2 = 1;
      mu = zeros(d, 1);
      for i = 1:d
        obs = find(~isnan(Y(:, i)));
        if length(obs)>0
          mu(i) = mean(Y(obs, i));
        else
          mu(i) = 0;
        end
      end
      numObs = sum(sum(~isnan(Y)));
      for i = 1:iters
        M = W'*W + sigma2*eye(q);
        invM = inv(M);
        exp_xxT = zeros(q);
        exp_x = zeros(N, q);
        for n = 1:N
          obs = find(~isnan(Y(n, :)));
          exp_x(n, :) = (invM*W(obs, :)'*(Y(n, obs)' - mu(obs)))';
        end
        exp_xxT = N*sigma2*invM + exp_x'*exp_x;
        s = zeros(d, q);
        s2 = 0;
        for n = 1:N
          obs = find(~isnan(Y(n, :)));
          subY = zeros(size(Y(n, :)))';
          subY(obs) = Y(n, obs)' - mu(obs);
          s = s + (subY)*exp_x(n, :);
          s2 = s2 + sum((Y(n, obs)' - mu(obs)).^2) - 2*exp_x(n, :)*W(obs, :)'*(Y(n, obs)' - mu(obs));
        end
        W = s*inv(exp_xxT);
        sigma2 = 1/(numObs)*(s2 + trace(exp_xxT*W'*W));
      end
      
      X = exp_x;
    end
def options = ppcaOptions(mappingType, latentPoints, mappingOptions):

    """Options for probabilistic PCA.
    
    Description:
    
    options = ppcaOptions returns the default options for
     probabilistic PCA.
     Returns:
      options - default options structure for probabilistic PCA.
        
        

    See also
    PPCACREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    options = [];
def [Y, Phi] = ppcaOut(model, X):

    """Output of an PPCA model.
    
    Description:
    
    Y = ppcaOut(model, X) gives the output of a density network for a
     given input.
     Returns:
      Y - the output.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
        

    See also
    PPCACREATE, MODELOUT


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    
    if nargin<2  
      % If we are just updating outer layer, basis functions are stored.
      Y = model.X*model.W + repmat(model.b, model.M, 1);
      if nargout>1
        Phi = model.X;
      end
    else
      Y = X*model.W + repmat(model.b, size(X, 1), 1);
      if nargout>1
        Phi = X;
      end
    end
def [mu, varsigma] = ppcaPosteriorMeanVar(model, X):

    """Mean and variances of the posterior at points given by X.
    
    Description:
    
    [mu, sigma] = ppcaPosteriorMeanVar(model, x) returns the posterior
     mean and variance for a given set of points.
     Returns:
      mu - the mean of the posterior distribution.
      sigma - the variances of the posterior distributions.
     Arguments:
      model - the model for which the posterior will be computed.
      x - the input positions for which the posterior will be computed.
        

    See also
    PPCACREATE


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    [mu, phi] = ppcaOut(model, X);
    %sqrt(det(phi*model.A*model.A'*phi))
    % Include magnification factors here ... need derivatives of outputs ...
    % modelOutputGradX ...
    varsigma = repmat(1/model.beta, size(X, 1), 1);

def model = rbfCreate(inputDim, outputDim, options):

    """Wrapper for NETLAB's rbf `net'.
    
    Description:
    def model = rbfCreate(inputDim, outputDim, options):
%
    """
        
    model = rbf(inputDim, options.hiddenDim, outputDim, options.activeFunc, options.outFunc);
    model.numParams = model.nwts;
    model.inputDim = inputDim;
    model.outputDim = outputDim;

def rbfDisplay(model, spacing):

    """Display an RBF network.
    
    Description:
    def rbfDisplay(model, spacing):
%
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Radial Basis Function network model:\n')
    fprintf(spacing);
    fprintf('  Input units: %d\n', model.inputDim);
    fprintf(spacing);
    fprintf('  Output units: %d\n', model.outputDim);
    fprintf(spacing);
    fprintf('  Hidden units: %d\n', model.nhidden);
    fprintf(spacing);
    fprintf('  Number of parameters: %d\n', model.numParams);
    fprintf(spacing);
    fprintf(['  Activation function: ' model.actfn '\n']);
    fprintf(['  Output function: ' model.outfn '\n']);
    

def model = rbfExpandParam(model, params):

    """Update rbf model with new vector of parameters.
    
    Description:
    
    model = rbfExpandParam(model, params) takes a vector of RBF
     weights and centres and places them in their respective positions
     in the RBF model. The function is a wrapper for the rbfunpak
     command.
     Returns:
      model - the model with the weights distributed in the correct
       places.
     Arguments:
      model - the model in which the weights are to be placed.
      params - a vector of the weights to be placed in the model.
        

    See also
    RBFUNPAK, RBFCREATE, RBFEXTRACTPARAM


    Copyright (c) 2006, 2007 Neil D. Lawrence
    
    """
        
    model = rbfunpak(model, params);
def [params, names] = rbfExtractParam(model):

    """Wrapper for NETLAB's rbfpak.
    
    Description:
    
    [params, names] = rbfExtractParam(model) returns a vector of all
     the weights and biases from a RBF network model. For single hidden
     layer models the function is a wrapper for the rbfpak command.
     Returns:
      params - vector of all the weights and biases returned by the
       model. The structure is governed by rbfpak.
      names - optional additional returned cell array of the names of
       the parameters.
     Arguments:
      model - the model from which we wish to extract the weights and
       biases.
        

    See also
    RBFPAK, RBFCREATE, RBFEXPANDPARAM, MODELEXTRACTPARAM


    Copyright (c) 2006, 2007, 2008 Neil D. Lawrence
    
    """
        
    
    params = rbfpak(model);
    if nargout > 1
      counter = 0;
      for j = 1:size(model.c, 2)
        for i = 1:size(model.c, 1)
          counter = counter + 1;
          names{counter} = ['Input centre ' num2str(i) '-' num2str(j)];
        end
      end
      for j = 1:size(model.wi, 2)
        counter = counter + 1;
        names{counter} = ['Hidden node width ' num2str(j)];
      end
      for j = 1:size(model.w2, 2)
        for i = 1:size(model.w2, 1)
          counter = counter + 1;
          names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
        end
      end
      for j = 1:size(model.b2, 2)
        counter = counter + 1;
        names{counter} = ['Output node bias ' num2str(j)];
      end
    end

def model = rbfOptimise(model, X, Y, display, iters):

    """Optimise RBF for given inputs and outputs.
    
    Description:
    def model = rbfOptimise(model, X, Y, display, iters):
%
    """
        
    if nargin < 4
      display = 1;
      if nargin < 5
        iters = 500;
      end
    end
    
    options = optOptions;
    options(14) = iters;
    options(1) = display;
    model = netopt(model, options, X, Y, 'scg');  
def options = rbfOptions

    """Default options for RBF network.
    
    Description:
    def options = rbfOptions
%
    """
        
    options.outFunc = 'linear';
    options.activeFunc = 'gaussian';
    options.hiddenDim = 20;

def [Y, G] = rbfOut(model, X):

    """Output of an RBF model.
    
    Description:
    
    Y = rbfOut(model, X) gives the output of a radial basis function
     model, the function is a wrapper for rbffwd.
     Returns:
      Y - the output.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
    
    [Y, G] = rbfOut(model, X) gives the output of a radial basis
     function model.
     Returns:
      Y - the output.
      G - the hidden layer activations.
     Arguments:
      model - the model for which the output is required.
      X - the input data for which the output is required.
        

    See also
    RBFFWD, RBF, MODELOUT


    Copyright (c) 2006, 2007, 2008 Neil D. Lawrence
    
    """
        
      if nargout > 1
        [Y, G] = rbffwd(model, X);
      else
        Y = rbffwd(model, X);
      end

def g = rbfOutputGrad(model, X):

    """Evaluate derivatives of rbf model outputs with respect to parameters.
    
    Description:
    
    g = rbfOutputGrad(model, X) evaluates the derivates of an RBF's
     outputs with respect to the parameters. Currently it simply wraps
     the NETLAB rbfderiv function.
     Returns:
      g - the gradient of the outputs of the RBF network with respect to
       each of the parameters. The size of the matrix is number of data x
       number of parameters x number of outputs of the model.
     Arguments:
      model - the model for which the derivatives are to be computed.
      X - the input data locations where the gradients are to be
       computed.
        

    See also
    RBFCREATE, RBFDERIV


    Copyright (c) 2006 Neil D. Lawrence
    
    """
        
    
    g = rbfderiv(model, X);

    """Evaluate derivatives of a RBF model's output with respect to inputs.
    
    Description:
    
    g = rbfOutputGradX(model, X) returns the derivatives of the
     outputs of an periodic radial basis function model with respect to
     the inputs to the model. Currently a wrapper for rbfjacob.
     Returns:
      g - the gradient of the output with respect to the inputs.
     Arguments:
      model - the model for which the derivatives will be computed.
      X - the locations at which the derivatives will be computed.
        

    See also
    RBFOUTPUTGRAD, MODELOUTPUTGRADX, RBFJACOB


    Copyright (c) 2008 Neil D. Lawrence
    
    """
        
    g = rbfjacob(model, X);

def model = rbfperiodicCreate(inputDim, outputDim, options):

    """Create a RBFPERIODIC model.
    
    Description:
        This model is a periodic, single input, model. The model is based on
        constructing an RBF network in the two dimensional space given by x_1 =
        sin(theta) and x_2 = cos(theta).  This leads to basis functions of
        the form.
        
        phi(theta) = exp(-2/sigma2*sin^2(0.5*(theta - m)))
        
        
    
    model = rbfperiodicCreate(inputDim, outputDim, options) creates a
     periodic radial basis function model structure given an options
     structure.
     Returns:
      model - the model structure with the default parameters placed in.
     Arguments:
      inputDim - the input dimension of the model.
      outputDim - the output dimension of the model.
      options - an options structure that determines the form of the
       model.
        

    See also
    RBFPERIODICKERNCREATE, RBFCREATE, RBFPERIODICOPTIONS, RBFPERIODICPARAMINIT, MODELCREATE


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    if inputDim>1
      error(['You may only create a one dimensional input periodic RBF ' ...
             'model.'])
    end
    model.inputDim = inputDim;
    model.outputDim = outputDim;
    
    model.widthTransform.type = optimiDefaultConstraint('positive');
    model.type = 'rbfperiodic';
    model.hiddenDim = options.hiddenDim;
    
    model.numParams = (inputDim+1)*options.hiddenDim + (options.hiddenDim + 1)*outputDim;
    model = rbfperiodicParamInit(model);
def rbfperiodicDisplay(model, spacing):

    """Display parameters of the RBFPERIODIC model.
    
    Description:
    
    rbfperiodicDisplay(model) displays the parameters of the periodic
     radial basis function model and the model type to the console.
     Arguments:
      model - the model to display.
    
    rbfperiodicDisplay(model, spacing)
     Arguments:
      model - the model to display.
      spacing - how many spaces to indent the display of the model by.
        

    See also
    RBFPERIODICCREATE, MODELDISPLAY


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    if nargin > 1
      spacing = repmat(32, 1, spacing);
    else
      spacing = [];
    end
    spacing = char(spacing);
    fprintf(spacing);
    fprintf('Periodic RBF model:\n')
    fprintf(spacing);
    fprintf('  Input dimensions: %d\n', model.inputDim);
    fprintf(spacing);
    fprintf('  Output dimensions: %d\n', model.outputDim);
    fprintf(spacing);
    fprintf('  Number of basis functions: %d\n', model.hiddenDim);
    fprintf(spacing);
    fprintf('  Number of parameters: %d\n', model.numParams);

def model = rbfperiodicExpandParam(model, params):

    """Create model structure from RBFPERIODIC model's parameters.
    
    Description:
    
    model = rbfperiodicExpandParam(model, param) returns a periodic
     radial basis function model structure filled with the parameters
     in the given vector. This is used as a helper function to enable
     parameters to be optimised in, for example, the NETLAB
     optimisation functions.
     Returns:
      model - model structure with the given parameters in the relevant
       locations.
     Arguments:
      model - the model structure in which the parameters are to be
       placed.
      param - vector of parameters which are to be placed in the model
       structure.
        

    See also
    RBFPERIODICCREATE, RBFPERIODICEXTRACTPARAM, MODELEXPANDPARAM


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        
    fhandle = str2func([model.widthTransform.type 'Transform']);
    startVal = 1;
    endVal = model.inputDim*model.hiddenDim; 
    model.thetaBar = reshape(params(startVal:endVal), model.inputDim, ...
                             model.hiddenDim);
    startVal = endVal+1;
    endVal = endVal + model.hiddenDim;
    model.sigma2 = reshape(fhandle(params(startVal:endVal), 'atox'), 1, model.hiddenDim);
    model.sigma2 = real(model.sigma2);
    startVal = endVal+1;
    endVal = endVal + model.hiddenDim*model.outputDim;
    model.weights = reshape(params(startVal:endVal), model.hiddenDim, ...
                            model.outputDim);
    startVal = endVal+1;
    endVal = endVal + model.outputDim;
    model.bias = reshape(params(startVal:endVal), 1, model.outputDim);

def [params, names] = rbfperiodicExtractParam(model):

    """Extract parameters from the RBFPERIODIC model structure.
    
    Description:
    
    param = rbfperiodicExtractParam(model) extracts parameters from
     the periodic radial basis function model structure into a vector
     of parameters for optimisation.
     Returns:
      param - vector of parameters extracted from the model.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
    
    [param, names] = rbfperiodicExtractParam(model) extracts
     parameters and parameter names from the periodic radial basis
     function model structure.
     Returns:
      param - vector of parameters extracted from the model.
      names - cell array of strings containing names for each parameter.
     Arguments:
      model - the model structure containing the parameters to be
       extracted.
        
        

    See also
    % SEEALSO RBFPERIODICCREATE, RBFPERIODICEXPANDPARAM, MODELEXTRACTPARAM, SCG, CONJGRAD


    Copyright (c) 2007 Neil D. Lawrence
    
    """
        fhandle = str2func([model.widthTransform.type 'Transform']);
    params = [model.thetaBar(:)' ...
    
