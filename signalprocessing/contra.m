%% function to contralateralise electrode data and collapse left/right sides

% data :    multi-dimensional matrix of data across all conditions
%           e.g. for 20 participants, 2 conditions, 2 sides, 64 electrodes,
%           size(data) = [20 2 2 64]
%           -- CANNOT be more than 20 dims (see code)
% sideVar : number of dimension containing side variable, e.g. 3
% elecVar : number of dimension containing electrode variable, e.g. 4
% keepSides : 1 or 0 (default 0), option to NOT collapse side variable

%


function data = contra(data,sideVar,elecVar,~)


    if nargin<4
        keepSides = 0;
    else
        keepSides = 1;
    end

    % check two levels to side variable
    if size(data,sideVar)~=2
        error('Not two levels to side variable')
    end
    
    % check at least 64 electrodes
    if size(data,elecVar)<64
        error('Not at least 64 electrodes')
    end
    
    % permute
    dims = size(data);
    vars = 1:length(dims);
    otherVars = find(~ismember(vars,[sideVar elecVar]));
    permuteOrder = [sideVar elecVar otherVars ];
    data = permute(data,permuteOrder);
    
    % contralateralise
    data(2,[  1:3    4:11  12:27  34:36  39:46   49:64 ],:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) ...
        = data(2,[ 34:36  39:46  49:64   1:3    4:11   12:27 ],:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
    
    % unpermute order
    [~,unpermuteOrder] = sort(permuteOrder);
    data = permute(data,unpermuteOrder);
    
    % collapse sides
    if ~keepSides
        data = collapse(data,sideVar);
    end
    
    
%     % EDITED 27/6/15
%     unpermuteOrder = [];
%     dimCount = 2;
%     for dim = 1:length(permuteOrder)
%         
%         if dim==sideVar
% %             dimCount = dimCount-1;
%             unpermuteOrder(dim) = 1;
%         elseif dim==elecVar
% %             dimCount = dimCount-1;
%             unpermuteOrder(dim) = 2;
%         else
%             dimCount = dimCount+1;
%             unpermuteOrder(dim) = dimCount;
%         end
%     end
%     
%     contraData = permute(contraData,unpermuteOrder);
%     
%     %collapse sides
%     contraData = collapse(contraData,sideVar);
    

end