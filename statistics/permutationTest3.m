%% permutation test with pixel- and cluster-wise corrections
% 
% adapted from Cohen ch.34, but separates positive and negative clusters
%
% inputs
%   cond1 & cond2   :   3D data sets, first two dimensions are x/y of map (e.g. frequency/time), with statistical 
%                       comparison performed across third dimension (e.g. elec x elec x trials, or elec x time x trials)
%
% versions: 
%   3. configured for cluster-based permutation, allow maxsum / maxsize cluster statistic
%   


% SHOULD AMALGAMATE THIS WITH CLUSTERPERM.M AND CLUSTERPERMPLOT.M


function [stat, h] = permutationTest3(cond1,cond2,varargin) % type,n_permutes,mcc_pval,voxel_pval)


% defaults
voxel_pval = .05;
mcc_pval = .025;
n_permutes = 1000;
correction_type = 'cluster';
cluster_statistic = 'maxsum';
plt = 0;
new_figure = 1;
savename = 0;
condition_labels = {'condition A','condition B'};
y_label = 1:size(cond1,2);
x_label = 1:size(cond1,1);
min_cluster_size = 1;
plot_data_type = 'zmap';


% input settings
for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'correction_type'
            correction_type = varargin{i+1};
        case 'n_permutes'
            n_permutes = varargin{i+1};
        case 'voxel_pval'
            voxel_pval = varargin{i+1};
        case 'mcc_pval'
            mcc_pval = varargin{i+1};
        case 'cluster_statistic'
            cluster_statistic = varargin{i+1};
        case 'plot'
            plt = varargin{i+1};
        case 'new_figure'
            new_figure = varargin{i+1};
        case 'condition_labels'
            condition_labels = varargin{i+1};
        case 'y_label'
            y_label = varargin{i+1};
        case 'x_label'
            x_label = varargin{i+1};
        case 'savename'
            savename = varargin{i+1};
        case 'min_cluster_size'
            min_cluster_size = varargin{i+1};
        case 'plot_data_type'
            plot_data_type = varargin{i+1};
        otherwise
            error('input not recognised')
    end
end



%% RUN PERMUTATIONS
if ~isstruct(cond1)

    % concatenate trials and map conditions
    combinedConds = cat(3,cond1,cond2);
    real_condition_mapping = [ ones(1,size(cond1,3)) -ones(1,size(cond2,3)) ];

    % compute actual t-test of difference (using unequal N and std)
    tnum   = squeeze(mean(combinedConds(:,:,real_condition_mapping==1),3) - mean(combinedConds(:,:,real_condition_mapping==-1),3)); 
    tdenom = sqrt( (std(combinedConds(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1) + (std(combinedConds(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) );
    real_t = tnum./tdenom;

    % initialize null hypothesis matrices
    permuted_tvals  = zeros(size(combinedConds,1),size(combinedConds,2),n_permutes);
    max_pixel_pvals = zeros(2,n_permutes);
    % max_clust_info  = zeros(n_permutes,1);

    % generate pixel-specific null hypothesis parameter distributions
    fprintf('permuting')
    max_clust_info = zeros(2,n_permutes);
    for permi = 1:n_permutes

        fprintf('.')
        fake_condition_mapping = sign(randn(size(combinedConds,3),1));

        % compute t-map of null hypothesis
        tnum   = squeeze(mean(combinedConds(:,:,fake_condition_mapping==1),3)-mean(combinedConds(:,:,fake_condition_mapping==-1),3));
        tdenom = sqrt( (std(combinedConds(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) + (std(combinedConds(:,:,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) );
        tmap   = tnum./tdenom;

        % save all permuted values
        permuted_tvals(:,:,permi) = tmap;

        % save maximum pixel values
        max_pixel_pvals(:,permi) = [ min(tmap(:)) max(tmap(:)) ];

        % threshold tmap to reveal clusters of significant pixels
        tmap(abs(tmap)<tinv(1-voxel_pval,size(combinedConds,3)-1))=0;

        % iterate separately for positive and negative clusters (departure from Mike Cohen code)
        cluster_valence{1} = max(0,tmap); % positive clusters
        cluster_valence{2} = min(0,tmap); % negative clusters
        for valence = 1:2
            
            % identify clusters
            clustinfo = bwconncomp(cluster_valence{valence}); 
            
            if strcmp(cluster_statistic, 'maxsize') 

                % get number of elements in largest supra-threshold cluster 
                max_clust_info(valence,permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps

            elseif strcmp(cluster_statistic, 'maxsum')

                % record largest sum of t-scores in cluster
                maxsum = zeros(1,length(clustinfo.PixelIdxList));
                for clust = 1:length(clustinfo.PixelIdxList)
                    maxsum(clust) = sum(abs(tmap(clustinfo.PixelIdxList{clust})));
                end
                if ~isempty(maxsum)
                    max_clust_info(valence,permi) = max(maxsum);
                end
                
            end
        end
        % ------------------------------------------------------------
        
    end
    fprintf('\n')

    
    % now compute Z-map
    zmap = (real_t - mean(permuted_tvals,3))./std(permuted_tvals,0,3);
        

    %% CORRECTION FOR MULTIPLE COMPARISONS

    % PIXEL-LEVEL CORRECTION
    if strcmp(correction_type,'pixel')

        % apply pixel-level corrected threshold
        lower_threshold = prctile(max_pixel_pvals(1,:),    mcc_pval*100/2);
        upper_threshold = prctile(max_pixel_pvals(2,:),100-mcc_pval*100/2);
        zmapthresh = zmap;
        zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=0; % pixel-level corrected threshold
        
        
    % CLUSTER-BASED CORRECTION
    elseif strcmp(correction_type,'cluster')

        % apply cluster-level corrected threshold
        zmapthresh = zmap;
        
        % uncorrected pixel-level threshold
        zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;

        % remove non-significant clusters from z-map
        if strcmp(cluster_statistic, 'maxsize') 
            
            % iterate separately for positive and negative clusters (departure from Mike Cohen code)
            cluster_valence{1} = max(0,zmapthresh); % positive clusters
            cluster_valence{2} = min(0,zmapthresh); % negative clusters
            for valence = 1:2

                % identify clusters
                clustinfo = bwconncomp(cluster_valence{valence}); 

                % find islands and remove those smaller than cluster size threshold
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                clust_threshold = prctile(max_clust_info(valence,:),100-mcc_pval*100);
            
                % identify clusters to remove
                whichclusters2remove = find(clust_info<clust_threshold);

                % remove clusters
                for i=1:length(whichclusters2remove)
                    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
                end
            end

        elseif strcmp(cluster_statistic, 'maxsum')

            % uncorrected pixel-level threshold
            tmapthresh = real_t;
            tmapthresh(abs(tmapthresh)<norminv(1-voxel_pval))=0; 
            
            % iterate separately for positive and negative clusters (departure from Mike Cohen code)
            cluster_valence{1} = max(0,tmapthresh); % positive clusters
            cluster_valence{2} = min(0,tmapthresh); % negative clusters
            cluster_pvals = cell(1,2);
            zmapthreshnew = zeros(size(zmapthresh));
            for valence = 1:2
                
                % identify clusters
                clustinfo = bwconncomp(cluster_valence{valence}); 
            
                % record sums of t-scores in clusters
                clust_info = zeros(1,length(clustinfo.PixelIdxList));
                clust_size = zeros(1,length(clustinfo.PixelIdxList));
                for clust = 1:length(clustinfo.PixelIdxList)
                    clust_info(clust) = sum(abs(tmapthresh(clustinfo.PixelIdxList{clust})));
                    clust_size(clust) = length(clustinfo.PixelIdxList{clust});
                end
                clust_threshold = prctile(max_clust_info(valence,:),100-mcc_pval*100);
                
                % keep significant clusters only 
                keepX = find(clust_info>clust_threshold & clust_size>=min_cluster_size);
                for i=1:length(keepX)
                    pixelIDX = clustinfo.PixelIdxList{keepX(i)};
                    zmapthreshnew(pixelIDX)=zmapthresh(pixelIDX);
                end
                zmapthresh = zmapthreshnew;
                
                % save significant p-vals
                for ii = 1:length(keepX)
                    all_scores = sort(max_clust_info(valence,:));
                    upper = find(clust_info(keepX(ii))<all_scores,1,'first');
                    lower = find(clust_info(keepX(ii))>all_scores,1,'last');
                    pval = 1-mean([upper,lower])/length(all_scores); 
                    cluster_pvals{valence} = [cluster_pvals{valence} pval];
                end
            end
        end
    end
    
    % save
    stat.zmap = zmap;
    stat.zmapthresh = zmapthresh;
    stat.voxel_pval = voxel_pval;
    stat.mcc_pval = mcc_pval;
    stat.n_permutes = n_permutes;
    stat.correction_type = correction_type;
    if strcmp(correction_type,'cluster')
        stat.cluster_statistic = cluster_statistic;
        stat.cluster_pvals = cluster_pvals;
    end
    if savename
        save(savename,'stat')
    end
    
else
    
    % load previously computed statistic
    stat = cond1;
    
end



%% PLOT
if plt

    if any(size(cond1)==1) 
        
        % line plot
        if new_figure
            h = figure('position',[200 100 800 400]); set(h,'Colormap',jet); hold on;
        else
            h=gcf; hold on;
        end
        plot(x_label,mean(cond1,3),'k-') % plot condition 1
        plot(x_label,mean(cond2,3),'k--')  % plot condition 2
        sigPixels = find(zmapthresh~=0); % plot significant differences
        plot(x_label(sigPixels), max(get(gca,'ylim'))*ones(size(sigPixels)), 'k.')
        legend(condition_labels)
        
    else

        % 2D plot
        if new_figure
            h = figure('position',[200 100 1200 400]);
            
            % plot condition 1
            subplot(1,3,1)
            imagesc(x_label,y_label,mean(cond1,3))
            lim = max(abs(get(gca,'clim')));
            set(gca,'clim',[-lim lim])
            colorbar
            title(condition_labels{1})

            % plot condition 2
            subplot(1,3,2)
            imagesc(x_label,y_label,mean(cond2,3))
            lim = max(abs(get(gca,'clim')));
            set(gca,'clim',[-lim lim])
            colorbar
            title(condition_labels{2})
            
            % set up axis for difference map
            subplot(1,3,3)
            
        end

        % plot difference map / zmap
        if strcmp(plot_data_type,'diff')
            h2 = imagesc(x_label,y_label,mean(cond1-cond2,3)); hold on % plot thresholded difference zmap into current axis
        elseif strcmp(plot_data_type,'zmap')
            h2 = imagesc(x_label,y_label,stat.zmap); hold on % plot thresholded difference zmap into current axis
%                 h2.AlphaData = 0.5 + 0.5*double(abs(stat.zmapthresh)>0);
        end
        lim = max(abs(get(gca,'clim')));
%         contour(x_label,y_label,double(logical(stat.zmapthresh)),1,'k')
%         contour(x_label,y_label,double(logical(stat.zmapthresh)),1,'w','linewidth',2)
        contour(x_label,y_label,double(logical(stat.zmapthresh)),1,'w')
        set(gca,'clim',[-lim lim])
        colorbar
        title([correction_type ' corrected difference zmap (p<' num2str(mcc_pval) ')'])

        % label clusters with p-values
        sign2 = {'+','-'};
        for ii = 1:2
            if ~isempty(stat.cluster_pvals{ii})
%                 text(x_label(round(2.1*length(x_label)/3)),y_label(round(length(y_label)/2)+ii),...
%                     ['p_' sign2{ii} '=' sprintf('%.3f',stat.cluster_pvals{ii}(1))]) % ARBITRARY LOCATION -- FIX
                text(x_label(round(2.1*length(x_label)/3)),y_label(round(length(y_label)/2)+ii),...
                    ['{\itp}=' sprintf('%.3f',stat.cluster_pvals{ii}(1))],'color','w') % ARBITRARY LOCATION -- FIX
            end
        end

    end
    
    % add title and save
    if savename
        saveasCS(h,savename,'png')
    end
    
    
        
%     h = figure('position',[0 0 1200 400]);
%     
%     % plot condition 1
%     subplot(1,3,1)
%     imagesc(x_label,y_label,mean(cond1,3))
%     lim = max(abs(get(gcf,'clim')));
%     set(gcf,'clim',[-lim lim])
%     colorbar
%     title(condition_labels{1})
%     
%     % plot condition 2
%     subplot(1,3,2)
%     imagesc(x_label,y_label,mean(cond2,3))
%     lim = max(abs(get(gcf,'clim')));
%     set(gcf,'clim',[-lim lim])
%     colorbar
%     title(condition_labels{2})
%     
%     % plot thresholded difference zmap
%     subplot(1,3,3)
%     h2 = imagesc(x_label,y_label,stat.zmap);
%     lim = max(abs(get(gcf,'clim')));
%     set(gcf,'clim',[-lim lim])
%     colorbar
%     a = 0.25 + 0.75*double(abs(stat.zmapthresh)>0);
%     h2.AlphaData = a;
%     title([correction_type ' corrected difference zmap'])
%     
%     % add title and save
%     if savename
%         saveasCS(h,savename,'png')
%     end
    
end    


