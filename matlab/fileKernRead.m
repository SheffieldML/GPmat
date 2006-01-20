function k = fileKernRead(kern, varargin)

% FILEKERNREAD Read kernel values from file or cache.

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

