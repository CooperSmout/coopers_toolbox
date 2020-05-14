function [im] = plotClusterPerm_chan_freq_time (clusterPermStat,fixedVariable,valueOfInterest,clusterAlpha,times2plot)

%% FUNCTION TO PLOT RESULTS OF SPATIOTEMPORAL CLUSTER-BASED PERMUTATION
% adapted from http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock
%
% only works for chan x freq x time cluster permutation
%


% variables
if nargin<4
    clusterAlpha = clusterPermStat.cfg.clusteralpha;
end

if strcmp(fixedVariable,'chan')
    
    data = squeeze(clusterPermStat.stat(valueOfInterest,:,:));
    figure
    plotTFR(data,clusterPermStat.time,clusterPermStat.freq)
    
    hold on
%     clust = squeeze(clusterPermStat.posclusterslabelmat(valueOfInterest,:,:));
    
    % mark clusters
    sign = {'pos','neg'};
    for ii = 1:2 % cycle through pos / neg clusters

        if isfield(clusterPermStat,[sign{ii} 'clusters']) && ~isempty(clusterPermStat.([sign{ii} 'clusters']))
            sig_clust = find([clusterPermStat.([sign{ii} 'clusters'])(:).prob] < clusterAlpha);
            dat = clusterPermStat.([sign{ii} 'clusterslabelmat']);
            dat = squeeze(dat(valueOfInterest,:,:));
            sig_clust_id = ismember(dat, sig_clust);

            if any(sig_clust_id(:))
                contour(sig_clust_id,'w')
            end
        end

    end
    
    
    
elseif strcmp(fixedVariable,'freq')
    
    % constrict data to selected frequency
    id = clusterPermStat.freq==valueOfInterest;
    clusterPermStat.diff = collapse(clusterPermStat.stat(:,id,:),2);
    clusterPermStat.posclusterslabelmat = collapse(clusterPermStat.posclusterslabelmat(:,id,:),2);
    clusterPermStat.negclusterslabelmat = collapse(clusterPermStat.negclusterslabelmat(:,id,:),2);
    
    % plot
    figure
    plotClusterPerm_chan_time (clusterPermStat,times2plot,clusterAlpha)
    suptitle([num2str(valueOfInterest) ' Hz'])
    
else
    error('must fix either channel or freq')
end


        