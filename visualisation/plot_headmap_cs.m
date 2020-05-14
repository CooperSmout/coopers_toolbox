function IM = plot_headmap( tit, head, cm, neg_limit, chan_locs, visible, dir_results, units ) 

h = figure('visible', visible);
if strcmp(neg_limit,'neg')
    limit = [ -max(abs(head)) max(abs(head))];
elseif strcmp(neg_limit,'pos')
    limit = [0 max(abs(head))];
else 
    limit = neg_limit;
end
    
topoplot_cs( head, chan_locs, 'maplimits', limit, 'colormap', cm, 'plotrad', .7, 'shading', 'interp', 'electrodes', 'labels' );
title( tit , 'interpreter', 'none' );
cbar_handle =  colorbar('peer', gca);
set( get( cbar_handle, 'title' ), 'String', units ); % , 'FontAngle', 'italic'

if strcmp(dir_results, 'no')
    disp('Image not saved')
else
    saveas(h, [ dir_results tit ], 'png'); 
    IM = imread( [ dir_results tit '.png'] );
end