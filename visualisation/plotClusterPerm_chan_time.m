function [im] = plotClusterPerm_chan_time (clusterPermStat,times2plot,filename,clusterAlpha,movieTiming,fig)


error('deprecated, use clusterPermChanTime instead')


%% FUNCTION TO PLOT RESULTS OF SPATIOTEMPORAL CLUSTER-BASED PERMUTATION
% adapted from http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock
%
% only works for chan x time cluster permutation
%
% edits:
% 05/02/16: 
    % changed to input actual times to plot (rather than timesteps)
    % added optional parameter plotAnySig to plot electrodes with ANY significant timepoint during plotted window (rather than ALL)
% 16/02/16:
    % added optional argument to change cluster alpha level (significance)
% 29/03/16 removed first input argument 'dif', cluster stat now requires additional field diff containing 
    % elec x time difference matrix (created when call ft_timelockstatisticsCS.m) 
% 18/04/16 added optional output argument to create frames of movie

   

% variables
if nargin<3
    filename=0;
elseif ~ischar(filename)
    error('filename for saving should be string')
end 
if nargin<4
    clusterAlpha = clusterPermStat.cfg.clusteralpha;
end
if nargin<5
    movieTiming=0;
end
if nargin<6
    fig=0;
end

% samples to plot
m = dsearchn(clusterPermStat.time',times2plot');

% subplot parameters
numSubplots = length(times2plot);
subs = [floor(sqrt(numSubplots)) ceil(sqrt(numSubplots))];
if subs(1)*subs(2)<numSubplots
    subs = [subs(2) subs(2)];
end

lim = max(abs(clusterPermStat.diff(:)));
maplimits = [-lim lim];

% figure for epoch plotting
if ~movieTiming
    im=figure('position',[0 0 800 600]); 
end

% movie parameters
movDim = [560 500];    % image dims
im3 = 255*ones(200,movDim(1),3,numSubplots-1); % snapshot window

% plot
for k = 1:numSubplots-1;

    % plot basic headmap
    dat = mean(clusterPermStat.diff(:,m(k):m(k+1)),2);
    
    
    if movieTiming 
        
        h4 = figure('position',[0 0 movDim]);
        
        currEndTime = times2plot(k+1);
        if ~mod(currEndTime,movieTiming) || k==1
            prevEndTime = currEndTime;
        end
%         title([ num2str(prevEndTime) ' ms'],'HorizontalAlignment','center')
        topoplot_cs(dat,'electrodes','off','maplimits',maplimits); 
%         axTopo = gca;
        set(gca,'fontsize',15)
        marksize = 12;
        
        % colorbar 
        cb=colorbar;
        set(cb,'position',[0.90    0.2489    0.0200    0.5368])
        text(.545,.38,'uV','fontsize',15) 
        
        
    else 
        
        subaxis(subs(1),subs(2),k); 
        title([ num2str(times2plot(k)) ' - ' num2str(times2plot(k+1)) ' ms'],'HorizontalAlignment','center')
        topoplot_cs(dat,'electrodes','off','maplimits',maplimits); 
        marksize = 5;
    
    end
    
    % plot significant clusters
    marks = {'*','o'};
    sign = {'pos','neg'};
    for ii = 1:2 % cycle through pos / neg clusters

        if isfield(clusterPermStat,[sign{ii} 'clusters']) && ~isempty(clusterPermStat.([sign{ii} 'clusters']))
            sig_clust = find([clusterPermStat.([sign{ii} 'clusters'])(:).prob] < clusterAlpha);
            sig_clust_id = ismember(clusterPermStat.([sign{ii} 'clusterslabelmat']), sig_clust);
            sig_clust_mark_any = any(sig_clust_id(:, m(k):m(k+1)), 2);
            sig_clust_mark_all = all(sig_clust_id(:, m(k):m(k+1)), 2);

            if any(sig_clust_mark_any)
                topoplot_cs(dat,'electrodes','off','maplimits',maplimits,'emarker2',{find(sig_clust_mark_any),marks{ii},'k',marksize,.1}); % plot significant 
            end
            if ~movieTiming && any(sig_clust_mark_all)
                topoplot_cs(dat,'electrodes','off','maplimits',maplimits,'emarker2',{find(sig_clust_mark_all),marks{ii},'k',8,1}); % plot significant 
            end
            maxCluster(ii) = clusterPermStat.([sign{ii} 'clusters'])(1).prob;
        
        else
            maxCluster(ii) = NaN;
        end

    end
    
    
    % draw time slice on figure
    if fig
        h5=open(fig); % open image
        y = ylim;
        if currEndTime>0
            line([currEndTime currEndTime],ylim,'color',[.5 .5 .5],'linewidth',4) % time slice
        end
        F2 = getframe(h5);
        axFig = gca;
        
        % capture image
        im2(:,:,:,k) = uint8(F2.cdata);
        close (h5)
        
    end
    
    
    % movie frames and timestamp
    if movieTiming
        
        % display time
        ax2 = axes('Position',[0 0 1 .9],'Visible','off'); % make new axis to enter time
        text(.5, .07, [ num2str(prevEndTime) ' ms'],'fontsize',20);
        
        % capture image
        F=getframe(h4);
        im(:,:,:,k) = uint8(F.cdata);
        close(h4)  
    end
    

end 

% join images
if fig
    im = cat(1, im2, im3, im);
end


% colorbar
if ~movieTiming
    subplot(subs(1),subs(2),numSubplots)
    topoplot_cs(zeros(1,64),'maplimits',maplimits,'electrodes','off')
    title(['p_c_r_i_t<' num2str(clusterPermStat.cfg.alpha)  ', p_c_l<' num2str(clusterAlpha)  ', cl_+=' sprintf('%.3f',maxCluster(1)) ', cl_-=' sprintf('%.3f',maxCluster(2))])
    colorbar
    
    % save image
    if filename
        saveasCS(im,filename,'png')
    end
    
end




end






%% OLD

%         % capture snapshot
%         if snapTimes
%             
%             % create figure with axes to plot into
%             pos = get(gcf,'position');
%             axPos = get(gca,'position');
%             axX = get(gca,'xlim');
% 
%             figSnap = figure('position',pos);
% %             axSnap = axes('position',axPos,'Visible','off','xlim',axX)
%             pos = get(gcf,'position');
%             set(gcf,'position',[pos(1:3) 200])
%         
%             % plot snapshot
%             if any(currEndTime==snapTimes)
%                 
%                 % turn snap on for subsequent plots
%                 snapOn(currEndTime==snapTimes)=1; 
%                 
%                 % new snap plot
%                 F=getframe(figSnap);
%                 snap = uint8(F.cdata);
% %                 snap = snap(1:4:end,1:4:end,:);
%                 im3(:,:,:,k) = snap;
%                 
%                 % line on plot
%                 for s = 1:sum(snapOn)
%                     axSnap = axes('position',[0 0 1 1],'Visible','off');
%                     copyobj(axTopo,figSnap)
%                     line([snapTimes(s) snapTimes(s)],[.95 1],'color',[.5 .5 .5],'linewidth',10) % time slice
%                 end
%                 
%             end
%         end
        