function blackPlot(xData, yData, SEM, colours, axes, saveName, varargin)

% data and SEM should be organised with matching conditions in row
% colours should be n x 3 matrix with RGB values for each condition
% legend is optional


figure('Position',[100 100 900 500]);
for cond = 1:size(yData,1)
    [h, semHandle] = boundedline(xData, yData(cond,:),SEM(cond,:)','alpha');
    set(h,'Color',colours(cond,:),'LineWidth',2.5)
    set(semHandle,'FaceColor',colours(cond,:))
    handles(cond) = h;
end
set(gcf,'Color','black','PaperPosition', [0 0 20 15])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'color','black','FontName','Arial','FontSize',20,'LineWidth',2,'xcolor','white','ycolor','white','box','off')
axis(axes)
% ylabel('Attentional Modulation (   dB power)','FontName','Arial','FontSize',25)
xlabel('Time (s)','FontName','Arial','FontSize',25)
line([0 length(xData)],[0 0],'linewidth',1,'linestyle','--','color','w')
if numel(varargin)
    legend(fliplr(handles),fliplr(varargin{1}),'FontName','Arial','FontSize',20,'textcolor','white')%,'position',[0.79 0.82 0.0688 0.0547])
end
% line([xData(64) xData(181)], [0.12 0.12],'color',[.3 .6 1])
% line([xData(64) xData(64)], [0.12 0.11],'color',[.3 .6 1])
% line([xData(181) xData(181)], [0.12 0.11],'color',[.3 .6 1])
% text(xData(120),.125,'*','color',[.3 .6 1],'fontsize',15)
% set(gca,'xticklabel',{'0' '' '1' '' '2' '' '3'})

saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

end



