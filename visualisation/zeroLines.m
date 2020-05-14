function zeroLines(ax,varargin)

% select axis
if nargin==0
    ax=gca;
end
axes(ax)
hold on

% if no input, plot lines at 0,0
if nargin<2
    varargin = {[0 0]};
end

if any(strcmp(varargin,':'))
    linestyle = ':';
    varargin  = varargin(~strcmp(varargin,':'));
else
    linestyle = '--';
end

% plot lines at specified locations
for zl = 1:length(varargin)
    x = varargin{zl}(1);
    y = varargin{zl}(2);
    xlim([min(x, min(get(ax,'xlim'))) max(x, max(get(ax,'xlim')))])
    ylim([min(y, min(get(ax,'ylim'))) max(y, max(get(ax,'ylim')))])
    plot([x x],get(ax,'ylim'),['k' linestyle])
    plot(get(ax,'xlim'),[y y],['k' linestyle])
end

