
function fig = plotTFR(data_tf,x,y,varargin)
% 
%   PLOTTFR
%   FIG = PLOTTFR(data_xy,x,y) plots data in xy format (e.g. time x freq) and labels axes
%   FIG = PLOTTFR(data_xy,x,y,'log10') plots log10 transformed data 
%   FIG = PLOTTFR(data_xy,x,y,'log2') plots log2 transformed data 
%   FIG = PLOTTFR(data_xy,x,y,'neg') forces clim to be +/-
% 
%   #########################################
%   #  Cooper Smout                         #
%   #  c.smout@uq.edu.au                    #
%   #  Queensland Brain Institute           #
%   #  University of Queensland, Australia  #
%   #########################################

fig = gcf;
data_yx = data_tf';


% LOG SCALING
if any(strcmp(varargin,'log10')) || any(strcmp(varargin,'log2'))  

    % set base of logarithm function
    if any(strcmp(varargin,'log10'))
        logBase = 10;
        logFunc = @log10;
    elseif any(strcmp(varargin,'log2'))
        logBase = 2;
        logFunc = @log2;
    end
    
    % scale positives and negatives separately 
    if any(data_yx(:)<=0) || any(strcmp(varargin,'neg'))
        
        % preallocate
        logData = zeros(size(data_yx));

        % arbitrary zero point
%         zeroPoint = min(logFunc(abs(data(:)))); % mn;
        zeroPoint = logFunc(10.^-6); % mn;
        if min(logFunc(abs(data_yx(:)))) < zeroPoint
            warning(['treating values below +/- 1e' num2str(zeroPoint) ' as zero'])
        end
        
        % transform positive values
        pos = data_yx>10.^zeroPoint;
        logData(pos) = logFunc(data_yx(pos)) - zeroPoint; % centre log data around zero

        % transform negative values
        neg = data_yx<10.^zeroPoint;
        logData(neg) = -logFunc(abs(data_yx(neg))) + zeroPoint ; % centre log data around zero
        
        % plot
%         contourf(time,freqs,logData,40,'linecolor','none')
        imagesc(x,1:length(y),logData)

        % modify colorbar
        cb = colorbar;
        logValues = get(cb,'ytick');
        posValues = logValues>0;
        negValues = logValues<0;
        logValues(posValues) = 10.^6 * logBase.^(logValues(posValues) + zeroPoint);   % convert log values to +mV
        logValues(negValues) = 10.^6 * -logBase.^(-logValues(negValues) + zeroPoint); % convert log values to -mV
        set(cb,'YTicklabel',logValues);
        title(cb,'log')
        
    else
        
        % plot
%         contourf(time,freqs,logFunc(data),40,'linecolor','none')
        imagesc(x,1:length(y),logFunc(data_yx))

        % modify colorbar
        cb = colorbar;
        logValues = get(cb,'ytick');
        logValues = logBase.^(logValues);   % convert log values to +mV
        set(cb,'YTicklabel',logValues);
        title(cb,'log')
        
    end
    
    
% NORMAL SCALING    
else 
    
    %     contourf(time,freqs,data,40,'linecolor','none')
    imagesc(x,1:length(y),data_yx);

    cb = colorbar;
%     title(cb,'mv')
    
    % limits
    if any(data_yx(:)<=0) || any(strcmp(varargin,'neg')) 
        set(gca,'clim',[-max(abs(get(gca,'clim'))) max(abs(get(gca,'clim')))]);
    else
        set(gca,'clim',[0 max(get(gca,'clim'))]);
    end
    
    
end

% plot settings
set(gca,'ydir','normal')
set(gca,'yticklabel',y)
set(gca,'tickdir','out')

% % contourf settings
% set(gca,'ytick',freqs)

% y labels
set(gca,'ytick',1:length(y))

% % x labels
% xdim = get(gca,'xlim');
% set(gca,'xtick',linspace(1,floor(xdim(2)),5))
% set(gca,'xtick',round(get(gca,'xtick')))
% set(gca,'xticklabel',x(get(gca,'xtick')))
