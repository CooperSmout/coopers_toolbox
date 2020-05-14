
function [h,ax] = figureCS(sz,varargin)

    % settings
    if nargin==0
        sz = [18 10];
    elseif ishandle(sz)
        h = sz; hold on;
        sz = h.PaperSize;
    elseif isnumeric(sz)
        h=figure; hold on;
    else
        error
    end
    
    % defaults
    h.Color = [1 1 1];
    h.PaperUnits = 'centimeters';
    h.Units = 'centimeters';
    h.PaperOrientation = 'landscape';
    h.PaperPosition = [5 5 sz];
    h.PaperSize = sz;
    h.PaperType = '<custom>';
    h.Position = [5 5 sz];
    h.Colormap = jet;
    h.PaperPositionMode = 'auto';
    
    % override defaults
    if nargin>1
        set(h,varargin{:})
    end
    
    % axes
    ax = gca;
    ax.FontName = 'Arial';
    ax.FontWeight = 'normal';
    ax.FontSize = 6;
    ax.TickDir = 'out';
    ax.LineWidth = 0.5; 
    ax.Box = 'off';
    
    % x-axis
    ax.XAxis.Color = [0 0 0];
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:.5:10;
    ax.XAxis.Label.String = 'XLabel';
    ax.XAxis.Label.FontSize = 8;
    ax.XAxis.Label.Color = [0 0 0];
    ax.XAxis.Limits = [0 10];
    
    % y-axis
    ax.YAxis.Color = [0 0 0];
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = 0:.5:10;
    ax.YAxis.Label.String = 'YLabel';
    ax.YAxis.Label.FontSize = 8;
    ax.XAxis.Label.Color = [0 0 0];
    ax.YAxis.Limits = [0 10];
    
%     % legend
%     l=legend('test');
%     l.FontName = 'Arial';
%     l.FontWeight = 'normal';
%     l.FontSize = 8;
%     l.Box = 'off';
    
    
    
    
    

end