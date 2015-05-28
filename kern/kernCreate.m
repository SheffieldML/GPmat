function kern = kernCreate(X, kernelType)

% KERNCREATE Initialise a kernel structure.
% FORMAT
% DESC creates a kernel matrix structure given an design matrix of
% DESC input points and a kernel type.
% ARG X : Input data values (from which kernel will later be computed).
% ARG type : Type of kernel to be created, some standard types are
% 'lin', 'rbf', 'white', 'bias' and 'rbfard'. If a cell of the form
% {'cmpnd', 'rbf', 'lin', 'white'} is used a compound kernel
% based on the sum of the individual kernels will be created. The
% 'cmpnd' element at the start of the sequence is
% optional. Furthermore, {'tensor', 'rbf', 'lin'} can be used to
% give a tensor product kernel, whose elements are the formed from
% the products of the two indvidual kernel's elements and
% {'multi', 'rbf', ...} can be used to create a block structured
% kernel for use with multiple outputs.
% Finally the form {'parametric', struct('opt1', {val1}), 'rbf'} can
% be used to pass options to other kernels.
% RETURN kern : The kernel structure.
%
% FORMAT
% DESC creates a kernel matrix structure given the dimensions of
% the design matrix and the kernel type.
% ARG dim : input dimension of the design matrix (i.e. number of features
% in the design matrix). 
% ARG type : Type of kernel to be created, as described above.
% RETURN kern : The kernel structure.
%
% SEEALSO : kernParamInit
%
% COPYRIGHT: Neil D. Lawrence, 2006
%
% MODIFICATIONS: Antti Honkela, 2009, Andreas Damianou, 2012, Carl
% Henrik Ek 2013
%
% KERN

if iscell(X)  
    for i = 1:length(X)
        dim(i) = size(X{i}, 2);
        %/~
        % If it's a cell, then we assume it is *not* containing the dimension of
        % the data.  --- Neil 23/2/2009
        %    if dim(i) == 1 & size(X{i}, 1) == 1;   % ???
        %       dim(i) = X{i};
        %     end
        %~/
    end
else
    % This is a bit of a hack to allow the creation of a kernel without
    % providing an input data matrix (which is sometimes useful). If the X
    % structure is a 1x1 it is assumed that it is the dimension of the
    % input data.
    dim = size(X, 2);
    if dim == 1 & size(X, 1) == 1;
        dim = X;
    end
end

if iscell(kernelType) && strcmp(kernelType{1}, 'parametric'),
    if ~isstruct(kernelType{2}) && strcmp(kernelType{2}, 'path'),
        kern = pathKernCreate(X,kernelType{3});
        kern = kernParamInit(kern);
        kern.Kstore = [];
        kern.diagK = [];

        return;
    end
    kern.options = kernelType{2};
    kernelType = kernelType{3};
end


if iscell(kernelType)
    kern.inputDimension = dim;
    switch kernelType{1}
      case 'multi'
        % multi output block based kernel.
        start = 2;
        kern.type = 'multi';
      case 'tensor'
        % tensor kernel type
        start = 2;
        kern.type = 'tensor';
      case 'cmpnd'
        % compound kernel type
        start = 2;
        kern.type = 'cmpnd';
      case 'translate'
        % translate kernel type
        start = 2;
        kern.type = 'translate';
      case 'velotrans'
        % velocity translate kernel type
        start = 2;
        kern.type = 'velotrans';
      case 'exp'
        % exponentiated kernel type
        start = 2;
        kern.type = 'exp';
      otherwise
        % compound kernel type
        start = 1;
        kern.type = 'cmpnd';
    end
    switch kern.type
      case 'multi'
        for i = start:length(kernelType)
            if iscell(X)
                kern.comp{i-start+1} = kernCreate(X{i-start+1}, kernelType{i});
                kern.diagBlockDim{i-start+1} = length(X{i-start+1});
            else
                kern.comp{i-start+1} = kernCreate(X, kernelType{i});
            end
            kern.comp{i-start+1}.index = [];
        end
      case {'tensor', 'cmpnd', 'translate', 'velotrans'}
        %%%---
        if strcmp(kernelType{1}, 'invcmpnd')
            kern.type = 'invcmpnd';
            start = 2;
            if iscell(X)
                
                % Every kernel has its own input space, length(dim) should be
                % equal to length(kernelType(start:end))
                if length(dim) ~= length(kernelType(2:end))
                    error('For the invcmpnd kernel a separate input domain must be given for every compound');
                end
                
                inds{1} = 1:dim(1);
                kern.comp{1} = kernCreate(X{1}, kernelType{2});
                kern.comp{1}.index = inds{1};
                for i = 2:length(kernelType)-1
                    lastInd = inds{i-1}(end);
                    inds{i} = lastInd+1:lastInd + dim(i) ;
                    kern.comp{i} = kernCreate(X{i}, kernelType{i+1});
                    kern.comp{i}.index = inds{i};
                end
            else
                for i = start:length(kernelType)
                    kern.comp{i-start+1} = kernCreate(X, kernelType{i});
                    kern.comp{i-start+1}.index = [];
                end
            end
        else
            %%%---
            for i = start:length(kernelType)
                kern.comp{i-start+1} = kernCreate(X, kernelType{i});
                kern.comp{i-start+1}.index = [];
            end
        end
      case 'exp'
        if start == length(kernelType)
            kern.argument = kernCreate(X, kernelType{start});
        else
            kern.argument = kernCreate(X, kernelType(start:end));
        end
    end
    kern = kernParamInit(kern);         % 
elseif isstruct(kernelType)
    % If a structure is passed, use it as the kernel.
    kern = kernelType;
else
    kern.type = kernelType;
    if iscell(X)
        kern.inputDimension = dim(i);
    else
        kern.inputDimension = dim;
    end
    kern = kernParamInit(kern);
end
kern.Kstore = [];
kern.diagK = [];
