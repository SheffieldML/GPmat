function mocapResultsCppBvh(fileName, bvhFileName, dataType, varargin)

% MOCAPRESULTSCPPBVH Load results from cpp file and visualise as a bvh format.
%
%	Description:
%
%	MOCAPRESULTSCPPBVH(FILENAME, BVHFILENAME, VISUALISATIONTYPE, ...)
%	loads results form a model file saved using C++ code and visualise.
%	 Arguments:
%	  FILENAME - the file saved by the C++ code.
%	  BVHFILENAME - the file containing the bvh data.
%	  VISUALISATIONTYPE - the root word associated with the
%	   visualisation.
%	  ... - optional additional arguments passed to the visualise
%	   commands.
%	
%
%	See also
%	FGPLVMVISUALISE, BVHREADFILE, FGPLVMREADFROMFILE


%	Copyright (c) 2005, 2006, 2007 Neil D. Lawrence


[model, labels] = fgplvmReadFromFile(fileName);
bvhStruct = bvhReadFile(bvhFileName);

maxInd = 0;
for i=1:length(bvhStruct.tree)
  for j= 1:size(bvhStruct.tree(i).rotInd, 2)
    if bvhStruct.tree(i).rotInd(j)>maxInd
      maxInd = bvhStruct.tree(i).rotInd(j);
    end
  end
  for j= 1:size(bvhStruct.tree(i).posInd, 2)
    if bvhStruct.tree(i).posInd(j)>maxInd
      maxInd = bvhStruct.tree(i).posInd(j);      
    end
  end
end
padding = maxInd - size(model.y, 2);
% Visualise the results
switch size(model.X, 2) 
 case 1
%  gplvmVisualise1D(model, [dataType 'Visualise'], [dataType 'Modify'], ...
%		   bvhStruct, padding, varargin{:});
  
 case 2
  fgplvmVisualise(model, labels, [dataType 'Visualise'], [dataType 'Modify'], ...
                 bvhStruct, padding, varargin{:});
  
 otherwise 
  error('No visualisation code for data of this latent dimension.');
end