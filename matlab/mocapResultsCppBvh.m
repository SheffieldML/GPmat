function mocapResultsCppBvh(fileName, bvhFileName, dataType, varargin)

% MOCAPRESULTSCPPBVH Load results from cpp file and visualise as a bvh format.

% MOCAP

[model, labels] = fgplvmReadFromFile(fileName);
bvhStruct = bvhReadFile(bvhFileName);

maxInd = 0;
for i=1:length(bvhStruct)
  for j= 1:size(bvhStruct(i).rotInd, 2)
    if bvhStruct(i).rotInd(j)>maxInd
      maxInd = bvhStruct(i).rotInd(j);
    end
  end
  for j= 1:size(bvhStruct(i).posInd, 2)
    if bvhStruct(i).posInd(j)>maxInd
      maxInd = bvhStruct(i).posInd(j);      
    end
  end
end
padding = maxInd - size(model.Y, 2);
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