


function oaxisCS(ax)

% new axes
x = get(gca,'xlim');
y = get(gca,'ylim');
axis off

% x axis
line(x,[0 0],'color','k') % x
xticks = get(gca,'xtick'); % major x ticks
xticks = -100:100:300; % TEMP REMOVE
xticks = xticks(xticks~=0);
xtickSize = max(abs(get(gca,'ylim')))/50;
for xtick = xticks
    line([xtick xtick],[-xtickSize xtickSize],'LineStyle','-','Color','k')
    text(xtick,-4*xtickSize,num2str(xtick),'HorizontalAlignment','center','fontsize',15)
end  

% minor x ticks
xMinortickSize = max(abs(get(gca,'ylim')))/100;
for xtick = [-50 50 150 250 350]
    line([xtick xtick],[-xMinortickSize xMinortickSize],'LineStyle','-','Color','k')
end  

% y axis
line([0 0],y,'color','k') % y
yticks = get(gca,'ytick');
yticks = -10:2:4; % TEMP REMOVE
yticks = yticks(yticks~=0);
ytickSize = max(abs(get(gca,'xlim')))/100;
for ytick = yticks
    line([-ytickSize ytickSize],[ytick ytick],'LineStyle','-','Color','k')
    text(-2*ytickSize,ytick,num2str(ytick),'HorizontalAlignment','right','fontsize',15)
end     

% minor y ticks
yMinortickSize = max(abs(get(gca,'xlim')))/200;
for ytick = -9:4
    line([-yMinortickSize yMinortickSize],[ytick ytick],'LineStyle','-','Color','k')
end 

% axis labels
text(x(2)+5*ytickSize,0,'ms','HorizontalAlignment','left','fontsize',15) % xlabel
text(0,y(2)+5*xtickSize,'uV','HorizontalAlignment','center','fontsize',15) % ylabel

