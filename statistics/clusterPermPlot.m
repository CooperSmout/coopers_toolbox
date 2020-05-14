

function [h,sigPosNeg] = clusterPermPlot(stat,varargin)

% CLUSTERPERMPLOT  Visualises significant differences computed with CLUSTERPERM.M
%
%   [h] = CLUSTERPERMPLOT(stat) plots data{1}, data{2}, data{1}-data{2}, and significant differences. 
%
%   default settings (see below) can be overridden with value-pair inputs, e.g.
%   [h] = CLUSTERPERMPLOT(stat,'plotRawData',0) plots data{1}-data{2} and significant differences. 
%   [h] = CLUSTERPERMPLOT(stat,'clusterMark','outline') plots significant differences as white outline. 
%   [h] = CLUSTERPERMPLOT(stat,'times2plot',0:50:500) plots headmaps at every 50ms up to 500ms (default is 100ms up to 300ms)
%
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)


%% SETTINGS

% defaults
cfg.alpha           = stat.cfg.alpha;
cfg.plotRawData     = 1;
cfg.plotDifference  = 1;
cfg.useCurrentAxis  = 0;
cfg.clusterMark     = 'alpha'; % or 'outline'
cfg.movieTiming     = 0;
cfg.clusterPval     = 1; % print p-values on plot: 0=none, 1=significant only, 2=all
cfg.conditions      = stat.conditions;
cfg.fixVar          = 'chan'; % ('chan_freq_time' only) which variable to fix for plotting ('chan' or 'time') 
cfg.colormap        = jet;
cfg.elecs           = 30;
cfg.times2plot      = 0:100:max(stat.time); % milliseconds
cfg.limit           = [];
cfg.clim            = [];
cfg.climDiff        = [];
cfg.ylim            = [];
cfg.withinSubjSem   = 0;
cfg.clusterCol      = 'k';

% overwrite defaults
for ii = 1:2:length(varargin)
    if ~isfield(cfg,varargin{ii})
        warning('unrecognized input')
    end
    cfg.(varargin{ii}) = varargin{ii+1};
end

% number of tested conditions
ncond = size(stat.data,3);

% checks
if strcmp(stat.dimord,'chan_time') && cfg.times2plot(end)>stat.time(end); error('invalid plot times'); end
if cfg.plotDifference && ncond>2
   warning('Cannot plot difference of >2 conditions, headmaps will show clusters only') 
   cfg.plotDifference = 0;
end
if cfg.useCurrentAxis; cfg.plotRawData = 0; end
if cfg.plotRawData && length(stat.conditions)==1
    warning('Plotting single condition only')
    cfg.plotRawData = 0;
end

% create figure?
if cfg.useCurrentAxis
    h=gcf;
else
    if strcmp(stat.dimord(1:4),'chan')
        h=figure('position',[100 100 800 800],'colormap',cfg.colormap);
    else
        h=figure('position',[100 100 400*(ncond+cfg.plotDifference) 300],'colormap',cfg.colormap);
    end
end
hold on

% plot settings

% find significant clusters
sigPosNeg = {false(size(stat.stat)),false(size(stat.stat))};
sign = {'pos','neg'};
% symbol = {'+','-'};
[pval,xy] = deal(cell(1,2));
for ii = 1:2 % cycle through pos / neg clusters
    if isfield(stat,[sign{ii} 'clusters']) && ~isempty(stat.([sign{ii} 'clusters']))
        
        % significant cluster p-values and locations
        sig_clust = find([stat.([sign{ii} 'clusters'])(:).prob] < cfg.alpha);
        if ~isempty(sig_clust)
            sigPosNeg{ii} = ismember(stat.([sign{ii} 'clusterslabelmat']), sig_clust); % mark significant clusters
        end
        
        % all cluster p-values and locations
        for cl = 1:length([stat.([sign{ii} 'clusters'])(:).prob])
            pval{ii}(cl) = stat.([sign{ii} 'clusters'])(cl).prob;
%             idx = find(any(stat.([sign{ii} 'clusterslabelmat'])==cl,1));
%             idy = find(any(stat.([sign{ii} 'clusterslabelmat'])==cl,2));
            idx = find(sum(stat.([sign{ii} 'clusterslabelmat'])==cl,2));
            idy = find(sum(stat.([sign{ii} 'clusterslabelmat'])==cl,1));
            xy{ii}(cl,:) = floor([median(idx) median(idy)]);
        end
    end
end

% plot data and mark significant clusters
switch stat.dimord

    case 'freq_time'
        
        % settings
        xorig = stat.timeOrig;
        xtested = stat.time;
        ylab = 'Frequency (Hz)';
        yorig = stat.freq;
        
        % plot settings
        [cfg,stat] = getLims(cfg,stat);
        
        % plot data
        if cfg.plotRawData 
            im = {};
            for ii = 1:ncond
                subplot(1,ncond+cfg.plotDifference,ii)
                im{ii} = imagesc(xorig,1:length(stat.freq),stat.data(:,:,ii));
                set(gca,'clim',cfg.clim)
                colorbar
                ylabel(ylab)
                xlabel('Time (s)')
                title(cfg.conditions{ii})
                set(gca,'ydir','normal')
                set(gca,'ytick',1:length(stat.freq))
                set(gca,'yticklabel',stat.freq)
            end
        end
        if cfg.plotDifference
            if ~cfg.useCurrentAxis
                subplot(1,3,3)
            end
            if size(stat.data,3)==1
                dif = stat.data;
            else
                dif = -diff(stat.data,[],3);
            end
            lim = max(abs(dif));
            im = {imagesc(xorig,1:length(stat.freq),dif)};
            set(gca,'clim',cfg.climDiff)
            colorbar
            ylabel(ylab)
            xlabel('Time (s)')
            title([cfg.conditions{1} ' - ' cfg.conditions{2}])
            set(gca,'ydir','normal')
            set(gca,'ytick',1:length(stat.freq))
            set(gca,'yticklabel',stat.freq)
        end
        
        % mark clusters
        if strcmp(cfg.clusterMark,'alpha')
            for ax = 1:length(im)
                im{ax}.AlphaData = 0.5*ones(length(yorig),length(xorig));
            end
        end
        for ii = 1:2
            if any(sigPosNeg{ii}(:))
                sig = false(length(yorig),length(xorig));
                XID = find(xorig==xtested(1)) : find(xorig==xtested(end));
                YID = find(yorig==ytested(1)) : find(yorig==ytested(end));
                sig(YID,XID) = sigPosNeg{ii};
                for ax = 1:length(im)
                    hold on
                    if strcmp(cfg.clusterMark,'alpha')
                        im{ax}.AlphaData(sig) = 1;
                    elseif strcmp(cfg.clusterMark,'outline')
                        axes(im{ax}.Parent); hold on
                        clim = get(gca,'clim');
                        contour(xorig,yorig,sig,1,cfg.clusterCol,'linewidth',2)
                        set(gca,'clim',clim);
                    end
                end
            else
                fprintf(['no ' sign{ii} ' clusters found\n'])
            end
            
            % print cluster p-values
            if cfg.clusterPval
                for cl = 1:length(pval{ii}) % label clusters with p-values
                    if pval{ii}(cl)<cfg.alpha % significant
                        text(xtested(xy{ii}(cl,1)),ytested(xy{ii}(cl,2)),[' {\itp}=' sprintf('%.3f',pval{ii}(cl))],'color',cfg.clusterCol,'HorizontalAlignment','center')
                    elseif cfg.clusterPval==2
                        text(xtested(xy{ii}(cl,1)),ytested(xy{ii}(cl,2)),[' {\itp}=' sprintf('%.3f',pval{ii}(cl))],'color',cfg.clusterCol,'HorizontalAlignment','center')
                    end
                end
            end
        end

        % mark test window
        if ~isequal(xorig,stat.time);
            hold on
            plot([stat.time(1) stat.time(1)],get(gca,'ylim'),'k:')
            plot([stat.time(end) stat.time(end)],get(gca,'ylim'),'k:')
            if strcmp(stat.dimord,'time_time')
                plot(get(gca,'xlim'),[stat.time(1) stat.time(1)],'k:')
                plot(get(gca,'xlim'),[stat.time(end) stat.time(end)],'k:')
            end
        end
        
        
    case {'time_time','circ_time'}
        
        % settings
        xorig = stat.timeOrig;
        xtested = stat.time;
        if strcmp(stat.dimord,'time_time')
            ylab = 'Time (s)';
            yorig = stat.timeOrig;
            ytested = stat.time;
        elseif strcmp(stat.dimord,'circ_time')
            ylab = 'Orientation (\circ)';
            ytested = stat.circ;
            if isfield(stat,'circOrig')
                yorig = stat.circOrig;
            else
                yorig = stat.circ;
            end
        end
        
        % plot settings
        [cfg,stat] = getLims(cfg,stat);
        
        % plot data
        if cfg.plotRawData 
            im = {};
            for ii = 1:ncond
                subplot(1,ncond+cfg.plotDifference,ii)
                im{ii} = imagesc(xorig,yorig,stat.data(:,:,ii));
                set(gca,'clim',cfg.clim)
                colorbar
                ylabel(ylab)
                xlabel('Time (s)')
                title(cfg.conditions{ii})
                set(gca,'ydir','normal')
            end
        end
        if cfg.plotDifference
            if ~cfg.useCurrentAxis
                subplot(1,3,3)
            end
            if size(stat.data,3)==1
                dif = stat.data;
            else
                dif = -diff(stat.data,[],3);
            end
            lim = max(abs(dif));
            im = {imagesc(xorig,yorig,dif)};
            set(gca,'clim',cfg.climDiff)
            colorbar
            ylabel(ylab)
            xlabel('Time (s)')
            title([cfg.conditions{1} ' - ' cfg.conditions{2}])
            set(gca,'ydir','normal')
        end
        
        % mark clusters
        if strcmp(cfg.clusterMark,'alpha')
            for ax = 1:length(im)
                im{ax}.AlphaData = 0.5*ones(length(yorig),length(xorig));
            end
        end
        for ii = 1:2
            if any(sigPosNeg{ii}(:))
                sig = false(length(yorig),length(xorig));
                XID = find(xorig==xtested(1)) : find(xorig==xtested(end));
                YID = find(yorig==ytested(1)) : find(yorig==ytested(end));
                sig(YID,XID) = sigPosNeg{ii};
                for ax = 1:length(im)
                    hold on
                    if strcmp(cfg.clusterMark,'alpha')
                        im{ax}.AlphaData(sig) = 1;
                    elseif strcmp(cfg.clusterMark,'outline')
                        axes(im{ax}.Parent); hold on
                        clim = get(gca,'clim');
                        contour(xorig,yorig,sig,1,cfg.clusterCol,'linewidth',1)
                        set(gca,'clim',clim);
                    end
                end
            else
                fprintf(['no ' sign{ii} ' clusters found\n'])
            end
            
            % print cluster p-values
            if cfg.clusterPval
                for cl = 1:length(pval{ii}) % label clusters with p-values
                    if pval{ii}(cl)<cfg.alpha % significant
                        text(xtested(xy{ii}(cl,1)),ytested(xy{ii}(cl,2)),[' {\itp}=' sprintf('%.3f',pval{ii}(cl))],'color',cfg.clusterCol,'HorizontalAlignment','center')
                    elseif cfg.clusterPval==2
                        text(xtested(xy{ii}(cl,1)),ytested(xy{ii}(cl,2)),[' {\itp}=' sprintf('%.3f',pval{ii}(cl))],'color',cfg.clusterCol,'HorizontalAlignment','center')
                    end
                end
            end
        end

        % mark test window
        if ~isequal(xorig,stat.time);
            hold on
            plot([stat.time(1) stat.time(1)],get(gca,'ylim'),'k:')
            plot([stat.time(end) stat.time(end)],get(gca,'ylim'),'k:')
            if strcmp(stat.dimord,'time_time')
                plot(get(gca,'xlim'),[stat.time(1) stat.time(1)],'k:')
                plot(get(gca,'xlim'),[stat.time(end) stat.time(end)],'k:')
            end
        end
        
        

    case 'time' % time series

        % plot time series
        hold on
        if cfg.withinSubjSem && ~isempty(stat.withinSubjSem)
            l=boundedline(stat.timeOrig,stat.data,repmat(stat.withinSubjSem,1,1,size(stat.data,2)),'alpha');
            fprintf('plotting within-subject sem')
        else
            sem(:,1,:) = stat.sem;
            l=boundedline(stat.timeOrig,stat.data,sem,'alpha');
            fprintf('plotting condition sem')
        end
        xlim([stat.timeOrig(1),stat.timeOrig(end)])
        if ~isempty(cfg.ylim)
            ylim(cfg.ylim)
        end
        legend(l,cfg.conditions,'location','eastoutside')
        
        % mark significant timepoints
        for ii = 1:2 % cycle through pos / neg clusters
            id = find(sigPosNeg{ii});
            plot(stat.time(id),(max(get(gca,'ylim')))*.95*ones(length(id)),'k.','MarkerSize',10)
%             for cl = 1:length(sig_pval{ii}) % label clusters with p-values
%                 text(stat.time(sig_xy{ii}(cl,1)), .9*max(get(gca,'ylim')), [' {\itp}=' sprintf('%.3f',sig_pval{ii}(cl))],'color','k')
%             end
        end
    
        % mark test window
        if ~isequal(stat.timeOrig,stat.time);
            plot([stat.time(1) stat.time(1)],get(gca,'ylim'),'k:')
            plot([stat.time(end) stat.time(end)],get(gca,'ylim'),'k:')
        end
        
        % label all clusters
        for ii = 1:2 % cycle through pos / neg clusters
            for cl = 1:length(pval{ii}) % label clusters with p-values
                text(stat.time(xy{ii}(cl,1)), (.95-.1*cl)*max(get(gca,'ylim')), [' {\itp}=' sprintf('%.3f',pval{ii}(cl))],'color',cfg.clusterCol,'HorizontalAlignment','center')
            end
        end
        
    
    case 'chan_time' 
        
        % load channel settings
        load biosemi64
        
        % samples to plot
        m = dsearchn(stat.timeOrig',cfg.times2plot');
        m2 = dsearchn(stat.time',cfg.times2plot');

        % subplot parameters
        numHeadmaps = length(cfg.times2plot)-1;
        subyx = [floor(sqrt(numHeadmaps)) ceil(sqrt(numHeadmaps))];
        if subyx(1)*subyx(2)<numHeadmaps
            subyx = [subyx(2) subyx(2)];
        end
        subyx(1) = subyx(1)+1;
        dif = -diff(stat.data,[],3);
        if isempty(cfg.clim)
            lim = max(abs(dif(:)));
            cfg.clim = [-lim lim];
        end
        mark = {'x','o'};

        % plot headmap data
        for k = 1:numHeadmaps;

            % plot difference headmap
%             subaxis(subxy(1),subxy(2),k); 
            subplot(subyx(1),subyx(2),k); 
            title([ num2str(cfg.times2plot(k)) ' - ' num2str(cfg.times2plot(k+1)) ' ms'],'HorizontalAlignment','center')
            if size(stat.data,3)>2
                dat = zeros(1,64); % plot zero headmap (cannot display difference of 3 conditions)
            else
                dat = mean(dif(:,m(k):m(k+1)),2);
            end
            topoplot(dat,biosemi64.chanlocs(1:64),'electrodes','off','maplimits',cfg.clim); hold on
%             topoplot_cs(dat,'electrodes','off','maplimits',cfg.clim); hold on

            % mark significant clusters
            marksize = 5;
            maxCluster = nan(1,2);
            for ii = 1:2 % positive and negative
                if isfield(stat,[sign{ii} 'clusters']) && ~isempty(stat.([sign{ii} 'clusters']))
                    sig_clust_mark_any = any(sigPosNeg{ii}(:, m2(k):m2(k+1)), 2);
                    if any(sig_clust_mark_any)
                        topoplot(dat,biosemi64.chanlocs(1:64),'electrodes','off','maplimits',cfg.clim,'emarker2',{find(sig_clust_mark_any),mark{ii},cfg.clusterCol,marksize,.1}); % plot significant 
                    end
                    sig_clust_mark_all = all(sigPosNeg{ii}(:, m2(k):m2(k+1)), 2);
                    if any(sig_clust_mark_all)
                        topoplot(dat,biosemi64.chanlocs(1:64),'electrodes','off','maplimits',cfg.clim,'emarker2',{find(sig_clust_mark_all),mark{ii},cfg.clusterCol,2*marksize,.1}); % plot significant 
                    end
                    maxCluster(ii) = stat.([sign{ii} 'clusters'])(1).prob;
                end
            end
        end

        % colorbar and p-values
        subplot(subyx(1),subyx(2),prod(subyx))
        topoplot(zeros(1,64),biosemi64.chanlocs(1:64),'maplimits',cfg.clim,'electrodes','off')
        title(['p_+=' sprintf('%.3f',maxCluster(1)) ', p_-=' sprintf('%.3f',maxCluster(2)) ' (\alpha_c_l=' num2str(cfg.alpha)  ', \alpha_p_x_l=' num2str(stat.cfg.clusteralpha) ')'])
        colorbar
        
        % plot time series
        subplot(subyx(1),3,[subyx(1)*3-2 subyx(1)*3-1]); hold on
        if cfg.withinSubjSem && ~isempty(stat.withinSubjSem)
            l=boundedline(stat.timeOrig,squeeze(mean(stat.data(cfg.elecs,:,:),1)),repmat(mean(stat.withinSubjSem(cfg.elecs,:),1),1,1,size(stat.data,3)),'alpha');
            fprintf('plotting within-subject sem')
        else
            sem(:,1,:) = mean(stat.sem(cfg.elecs,:,:),1);
            l=boundedline(stat.timeOrig,squeeze(mean(stat.data(cfg.elecs,:,:),1)),sem,'alpha');
            fprintf('plotting condition sem')
        end
        legend(l,cfg.conditions,'location','eastoutside')
        axis tight
        plot(get(gca,'xlim'),[0 0],'k')
        plot([0 0],get(gca,'ylim'),'k')
        if ~isempty(cfg.ylim)
            ylim(cfg.ylim)
        end
        
        % mark pos / neg clusters
        for ii = 1:2 
            id = find(any(sigPosNeg{ii}(cfg.elecs,:),1));
            plot(stat.time(id),(max(get(gca,'ylim')))*ones(length(id)),'k.','MarkerSize',10)
        end
        title(biosemi64.elecs(cfg.elecs))
        
    case 'chan_freq_time' 

        error('Not coded for this dimord (for starter code, see plotTFR.m and /Users/uqcsmout/Projects/rp2/scripts/rp2_normtf.m')
        
        
    otherwise
        
        error('Not coded for this dimord (see code below for superceded chan-freq code')
    
    
end



%% ########################## functions
function [cfg,stat] = getLims(cfg,stat)

% check for and compute regular limits
if isempty(cfg.clim)
    cfg.clim = [min(stat.data(:)) max(stat.data(:))];
end

% check for and compute difference limits
if isempty(cfg.climDiff)
    diffs = diff(stat.data,[],3);
    lim = max(abs(diffs(:)));
    cfg.climDiff = [-lim lim];
end


















%% SUPERCEDED CODE USED TO MAKE MOVIES
% 
% % movie parameters
% movDim = [560 500];    % image dims
% im3 = 255*ones(200,movDim(1),3,numSubplots-1); % snapshot window
% 
% h4 = figure('position',[0 0 movDim]);
% currEndTime = cfg.times2plot(k+1);
% if ~mod(currEndTime,cfg.movieTiming) || k==1
%     prevEndTime = currEndTime;
% end
% %         title([ num2str(prevEndTime) ' ms'],'HorizontalAlignment','center')
% topoplot_cs(dat,'electrodes','off','maplimits',cfg.limit); 
% %         axTopo = gca;
% set(gca,'fontsize',15)
% marksize = 12;
% 
% % colorbar 
% cb=colorbar;
% set(cb,'position',[0.90    0.2489    0.0200    0.5368])
% text(.545,.38,'uV','fontsize',15) 
% 
% 
% % draw time slice on figure
% if fig
%     h5=open(fig); % open image
%     y = ylim;
%     if currEndTime>0
%         line([currEndTime currEndTime],ylim,'color',[.5 .5 .5],'linewidth',4) % time slice
%     end
%     F2 = getframe(h5);
%     axFig = gca;
%     im2(:,:,:,k) = uint8(F2.cdata); % capture image
%     close (h5)
% 
% end
% 
% % movie frames and timestamp
% if cfg.movieTiming
%     ax2 = axes('Position',[0 0 1 .9],'Visible','off'); % make new axis to enter time
%     text(.5, .07, [ num2str(prevEndTime) ' ms'],'fontsize',20); % display time
%     F=getframe(h4); % capture image
%     im(:,:,:,k) = uint8(F.cdata);
%     close(h4)  
% end
% 
% % join images
% if fig
%     im = cat(1, im2, im3, im);
% end



%% SUPERCEDED CODE TO PLOT CHAN-FREQ STATISTICS
% 
% if strcmp(fixedVar,'chan') % COLLAPSE ELECTRODES AND PLOT TIME X FREQ
% 
%     % collapse electrodes and plot
%     dat_xy = collapse(stat.diff(fixedVals,:,:),1)';
%     plotTFR(dat_xy,stat.timeOrig,stat.freqOrig)
% 
%     % work out significant freq x time clusters across plotted electrodes
%     cols = {'k','w'};
%     sign = {'pos','neg'};
%     for ii = 1:2 % cycle through pos / neg clusters
% 
%         if isfield(stat,[sign{ii} 'clusters']) && ~isempty(stat.([sign{ii} 'clusters']))
%             sig_clust = find([stat.([sign{ii} 'clusters'])(:).prob] < clusterAlpha);
%             sig_clust_id = ismember(stat.([sign{ii} 'clusterslabelmat']), sig_clust);
%             sig_clust_xy{ii} = squeeze( any(sig_clust_id(fixedVals,:,:), 1))'; 
% %                 sig_clust_xy{ii} = squeeze( all(sig_clust_id(fixedVals,:,:), 1))';
%             maxCluster(ii) = stat.([sign{ii} 'clusters'])(1).prob;
% 
%         else
%             sig_clust_xy{ii} = false(size(stat.prob,3),size(stat.prob,2));
%             maxCluster(ii) = NaN;
%         end
% 
% 
% %             % PLOT CONTOURS - NOT WORKING
% %             if any(sig_clust_any_TF{ii}(:))
% %             
% %                 % work out contours
% %                 C = contourc(double(sig_clust_any_TF{ii}),[1 1]);
% %                 clusterCount = 0;
% %                 clusterID = 0;
% %                 for cl = 1:length(C)
% %                     if clusterCount==0
% %                         clusterCount = C(2,cl);
% %                         clusterID = clusterID +1;
% %                         clusters{clusterID} = [];
% %                         clusterVal(clusterID) = C(1,cl);
% %                     else
% %                         clusters{clusterID} = [clusters{clusterID} C(:,cl)];
% %                         clusterCount = clusterCount-1;
% %                     end
% % 
% %                 end
% % 
% %                 % plot contours over TFR
% %                 hold on
% %                 for cont = 1:length(clusters)
% %                     f = fill(dat1.time(clusters{cont}(2,:)),dat1.freq(clusters{cont}(1,:)),cols{ii},'EdgeColor',cols{ii});
% %                     set(f,'FaceAlpha',.2)
% %                 end
% % 
% %             end
%     end
% 
% %         % threshold data and plot
% %         sig_clust_mark_any_both = sig_clust_any_TF{1} | sig_clust_any_TF{2};
% %         if any(sig_clust_mark_any_both(:))
% %             dat(~sig_clust_mark_any_both) = 0;
% %             subplot(1,2,2)
% %             plotTFR(dat,dat1.time,dat1.freq)
% % %             imagesc(dat1.time,dat1.freq,dat')
% %             title('thresholded')
% %         end
% 
% 
% elseif strcmp(fixedVar,'freq') % COLLAPSE FREQUENCIES AND PLOT ELEC X TIME
% 
%     % collapse frequencies and plot
%     dat_xy = collapse(stat.diff(:,fixedVals,:),2)';
% 
%     % rearrange biosemi electrode order to something more logical
%     order = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 33 37 38 47 48 ...
%         32 31 30 29 28 34 35 36 39 40 41 42 43 44 45 46 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64];
%     dat_xy = dat_xy(:,order);
%     plotTFR(dat_xy,stat.timeOrig,1:64)
%     set(gca,'yticklabel',biosemi64.elecs(order));
% 
%     % work out significant freq x time clusters across plotted electrodes
%     cols = {'k','w'};
%     sign = {'pos','neg'};
%     for ii = 1:2 % cycle through pos / neg clusters
% 
%         if isfield(stat,[sign{ii} 'clusters']) && ~isempty(stat.([sign{ii} 'clusters']))
%             sig_clust = find([stat.([sign{ii} 'clusters'])(:).prob] < clusterAlpha);
%             sig_clust_id = ismember(stat.([sign{ii} 'clusterslabelmat']), sig_clust);
%             sig_clust_xy{ii} = squeeze( any(sig_clust_id(:,fixedVals,:), 2))'; 
% %                 sig_clust_xy{ii} = squeeze( all(sig_clust_id(:,fixedVals,:), 2))'; 
%             maxCluster(ii) = stat.([sign{ii} 'clusters'])(1).prob;
% 
%         else
%             sig_clust_xy{ii} = false(size(stat.prob,3),size(stat.prob,1));
%             maxCluster(ii) = NaN;
%         end
%     end
% 
% end
% 
% % mark significant clusters
% nonsigAlpha = 0.3;
% if isfield(stat,'latencyIDs')
%     latencyIDs = stat.latencyIDs;
% else
%     latencyIDs = [1 size(stat.prob,3)];
% end
% if isfield(stat,'frequencyIDs')
%     frequencyIDs = stat.frequencyIDs;
% else
%     frequencyIDs = [1 size(stat.prob,2)];
% end
% sig_clust_both_yx = zeros(size(dat_xy'));
% sig_clust_both_yx(frequencyIDs(1):frequencyIDs(2),latencyIDs(1):latencyIDs(2)) = double(sig_clust_xy{1} | sig_clust_xy{2})';
% 
% % draw lines to mark tested time-frequency points
% x1 = stat.timeOrig(stat.latencyIDs(1));
% x2 = stat.timeOrig(stat.latencyIDs(2));
% line([x1 x1],frequencyIDs,'color','w','LineWidth',2) % cluster perm start
% line([x2 x2],frequencyIDs,'color','w','LineWidth',2) % cluster perm end
% line([x1 x2],[frequencyIDs(2) frequencyIDs(2)],'color','w','LineWidth',2) % cluster perm top
% line([x1 x2],[frequencyIDs(1) frequencyIDs(1)],'color','w','LineWidth',2) % cluster perm bottom
% 
% alpha(double(sig_clust_both_yx + nonsigAlpha))
% 


