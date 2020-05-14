%% Function to rearrange electrodes and plot connectivity matrix
% Cooper Smout: c.smout@uq.edu.au
% Queensland Brain Institute, University of Queensland, Australia
%
% dat: connectivity matrix (elecs x elecs) with top right (obsolete) half of matrix as nans


function data = plotConnectivity(data,limits,notContra)

% set defaults
if nargin<3
   notContra = 0; 
   
    if nargin<2
        limits=[0 1];
    end
end

% check data
numElecs=size(data,1);
if numElecs~=size(data,2)
    error('matrix not square')
elseif ndims(data)>2
    error('too many dimensions in data')
end
% for E1 = 1:numElecs
%     for E2 = E1+1:numElecs
%         if abs(data(E1,E2))>0
%             flag = 1;
%         end
%     end
% end
% if exist('flag')
%     warning('values present in top half of data')
% end

% rearrange biosemi electrode order to something more logical
order    = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 33 37 38 47 48 ...
                32 31 30 29 28 34 35 36 39 40 41 42 43 44 45 46 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64];
if notContra
    label  	= {'LF'   'LC'    'LP'    'F'     'C'     'P'     'RF'    'RC'    'RP'};
else
    label   = {'IF'   'IC'    'IP'    'F'     'C'     'P'     'CF'    'CC'    'CP'};
end
labelTick   = [4.0    13.5    23.5    29.0    32.0    35.5    41.0    50.5    60.5];
boundaries  = [    7.5    19.5    27.5    30.5    33.5    37.5    44.5    56.5];

% permute electrodes
data = data(order,order);

% move flipped values back into matrix
data = correctConnectivity(data);

% plot
imagesc(data,limits);
set(gca,'xAxisLocation','top')
set(gca,'ticklength',[0 0],'xtick',labelTick,'xticklabel',label,'ytick',labelTick,'yticklabel',label,'ydir','normal')
axis square
colorbar

% electrode line dividers
for L = 1:length(boundaries)
    line([boundaries(L) boundaries(L)],[0 numElecs+1],'color','k','linewidth',1); % vertical lines
    line([0 numElecs+1],[boundaries(L) boundaries(L)],'color','k','linewidth',1); % horizontal lines
end




