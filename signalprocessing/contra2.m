%% function to contralateralise electrode data and collapse left/right sides

% data :    multi-dimensional matrix of data across all conditions
%           e.g. for 20 participants, 2 conditions, 2 sides, 64 electrodes,
%           size(data) = [20 2 2 64]
%           -- CANNOT be more than 20 dims (see code)
% sideDim : number of dimension containing side variable, e.g. 3
% elecDims : number of dimension containing electrode variable, e.g. 4 (default = 1)
%           -- can also contain two dimensions if contralateralising a connectivity matrix
%
% if want to keep side variable (or doesn't exist), only input data and electrode dimension)
% if only input data, will assume that electrodes are on first dimension, and if side exists it is on second dimension


function data = contra2(data,elecDims,sideDim,keepSides)

    % keep side dimension?
    if nargin<4
        keepSides = 0;
    end
    
    % if only input data and electrode dimension, assumes data is single-sided
    if nargin<3
        sideDim = ndims(data)+1; % suitably large number that it will just pick single dimension
    end

    % check side dimension
    if size(data,sideDim)==1
        fprintf('Data inputted as single-sided...')
        side2flip = 1;
    elseif size(data,sideDim)>2
        error('More than two levels to side dimension')
    else
        fprintf('Data is two-sided, will flip second side...')
        side2flip = 2; % flip right side
    end
    
    % electrode dimensions
    if nargin<2
        elecDims = 1; % assumes electrodes on first dimension if dimension not given
        warning('Electrode dimension not entered, flipping values in first dimension...')
    elseif length(elecDims)>2
        error('too many electrode dims')
    elseif length(elecDims)==2
        fprintf('Flipping connectivity matrix...')
    elseif length(elecDims)==1
        fprintf('Flipping electrode headmap...')
    end
    
    
    % contralateralise each electrode dimension
    for elecVar = elecDims
        
         % check electrodes
        if size(data,elecVar) ~= 64 
            error('Not 64 electrodes!')
        end
        
        
        % permute data
        dims = size(data);
        vars = 1:length(dims);
        otherVars = find(~ismember(vars,[sideDim elecVar]));
        permuteOrder = [elecVar sideDim otherVars ];
        data = permute(data,permuteOrder);

        % contralateralise
        data([  1:3    4:11  12:27  34:36  39:46   49:64 ],side2flip,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) ...
            = data([ 34:36  39:46  49:64   1:3    4:11   12:27 ],side2flip,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:); 

        % unpermute order
        [~,unpermuteOrder] = sort(permuteOrder);
        data = permute(data,unpermuteOrder);

        % collapse side dimension? 
        if length(elecDims)==1 && ~keepSides % on first run through non-correlation matrix
            data = collapse(data,sideDim); % collapse side dimension
            
        % correct connectivity matrix?
        elseif length(elecDims)==2 && elecVar==elecDims(2) % on second run through correlation matrix
            
            % fix correlation matrix
            data = correctConnectivity(data); 
            
            if ~keepSides % collapse side dimension
                data = collapse(data,sideDim);
            end
            
            
        end
        
        
    end
    
  fprintf('Done!\n')

end






%% TEST THIS CODE IS DOING WHAT IT SHOULD (COPY TO NEW WINDOW!!)
% 
% % connectivity matrix
% orig = [];
% for E1 = 1:64
%     for E2 = 1:E1
%         if E1==E2
% %             orig(E1,E2)=1;
%             orig(E1,E2,:)=[1 1 1 1 1];
%         else
% %             orig(E1,E2)= E1*.015+E2*.0015;
%             orig(E1,E2,:)= (E1*.015+E2*.0015) * [1 1 1 1 1];
%         end
%     end
% end
% flipped = contra2(orig,[1 2]); % contralateralise
% figure; imagesc(orig(:,:,1),[0 1])
% figure; imagesc(flipped(:,:,1),[0 1])
% 
% 
% % headmap
% head = rand(64,1).*(1:64)';
% figure; topoplot_cs(head)
% figure; topoplot_cs(contra2(head))

