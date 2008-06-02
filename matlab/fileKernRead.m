function k = fileKernRead(kern, varargin)

% FILEKERNREAD Read kernel values from file or cache.
% FORMAT
% DESC returns values of a kernel which is stored on disk. On the
% first call the kernels are read into memory. On further calls,
% the cache is accessed.
% ARG kern : the kernel structure for which file is to be returned.
% ARG index : indices of the kernel which are to be returned.
% RETURN k : elements of the kernel that were requested.
%
% FORMAT 
% DESC reads in kernel files with row and column indices specified.
% ARG kern : the kernel structure specifying the file to be loaded
% in.
% ARG index1 : row indices of the kernel which are to be returned.
% ARG index2 : column indices of the kernel which are to be returned.
% RETURN k : elements of the kernel that were requested.
%
% FORMAT 
% DESC reads in the diagonal of a kernel file with indices specified.
% ARG kern : the kernel structure specifying the file to be loaded
% in.
% ARG index : indices of the diagonal which are to be returned.
% ARG flag : set to 'diag' to return diagonal elements.
% RETURN k : elements of the kernel that were requested.
%
% SEEALSO : fileKernParamInit

% KERN

% cache for the kernel.
persistent K
persistent diagK
persistent index
persistent fileNames

loadFile = 0;
if isempty(K)
  if ~isfield(kern, 'fileName')
    error('Missing fileName field in file kernel.')
  end
  loadFile = 1;
  index = 1;
end

curIndex = -1;
for i = 1:length(fileNames)
  if strcmp(kern.fileName, fileNames{i})
    curIndex = i;
    break
  end
end
if curIndex == -1
  loadFile = 1;
end

if loadFile
  fprintf('Loading kernel ... %s\n', kern.fileName); 
  K{index + 1} = single(load(kern.fileName));
  fprintf('Loaded kernel %s into memory.\n', kern.fileName); 
  index = index+1;
  diagK{index} = diag(K{index});
  fileNames{index} = kern.fileName;
  curIndex = index;
end

if isstr(varargin{end})
  if strcmp(varargin{end}, 'diag')
    k = double(diagK{curIndex}(varargin{1}));
    return
  end
end

if length(varargin)>1
  k = double(K{curIndex}(varargin{1}, varargin{2}));
else
  k = double(K{curIndex}(varargin{1}, varargin{1}));
end

