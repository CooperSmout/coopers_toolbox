function [pnts,hs] = clusterPermMark(stat,elec,txt,ysig,color,linestyle)

% CLUSTERPERMMARK  Plots significant timepoints to figure created with with CLUSTERPERMPLOT.M
%
%   edits:
%       09 Oct 2018: Changed format for ysig input to be a percentage of ylim off lower axis
%       02 Apr 2019: Included option to not plot if text==nan --- NOT A E
%
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)


% defaults
if nargin<3
    txt=0;
end
if nargin<4 || isempty(ysig)
%     yl = get(gca,'ylim');
%     ysig = yl(1)+diff(yl)/20;
    ysig = .05;
end
if nargin<5 || isempty(color)
    color=[0 0 0];
end
if nargin<6 || isempty(linestyle)
    linestyle = '-';
end
yl = get(gca,'ylim');
ysig = ysig*diff(yl) + min(yl);
ytxt = ysig+diff(yl)/20;

hold on
signs = {'pos','neg'};
pnts = false(1,length(stat.prob));
hs=[];
for val = 1:2
    fld = [signs{val} 'clusters'];
    if isfield(stat,fld) && ~isempty(stat.(fld))
        probs = [stat.([signs{val} 'clusters']).prob];
        id = find(probs<=.05);
        for ii = 1:length(id)
            if size(stat.prob,2)>1
                sig = any(stat.([signs{val} 'clusterslabelmat'])(elec,:) == id(ii),1);
            else
                sig = stat.([signs{val} 'clusterslabelmat']) == id(ii);
            end
%             pnts = [pnts find(sig)];
            pnts(sig) = true;
            if any(sig) && ~isnan(txt)
                cls = bwconncomp(sig);
                for cl = 1:length(cls.PixelIdxList)
%                     x1 = 1000*min(stat.time(cls.PixelIdxList{cl}));
%                     x2 = 1000*max(stat.time(cls.PixelIdxList{cl}));
                    x1 = min(stat.time(cls.PixelIdxList{cl}));
                    x2 = max(stat.time(cls.PixelIdxList{cl}));
%                     hs(end+1) = plot([x1 x2],[ysig ysig],'k-','linewidth',2,'color',color,'linestyle',linestyle);
                    hs(end+1) = plot([x1 x2],[ysig ysig],'k-','linewidth',2.5,'color',color,'linestyle',linestyle);
                end
                if txt
                    txt(mean(1000*stat.time(sig)),ytxt,['{\itp} = ' num2str(probs(id(ii)),'%.3f')],...
                    'horizontalalignment','center','FontName','Arial','FontWeight','normal','FontSize',9);
                end
            end
        end
    end
end
