function blackPlotWaveletCentred(xData, yData, SEM, colours, axes, saveName, varargin)

% data and SEM should be organised with matching conditions in row
% colours should be n x 3 matrix with RGB values for each condition
% legend is optional

if numel(varargin)
    lgnd = varargin;
else
    lgnd = 0;
end

figure('Position',[100 100 900 500]);
for cond = 1:size(yData,1)
    [h, semHandle] = boundedline(xData, yData(cond,:),SEM(cond,:)','alpha');
    set(h,'Color',colours(cond,:),'LineWidth',2.5)
    set(semHandle,'FaceColor',colours(cond,:))
    handles(cond) = h;
end
set(gcf,'Color','black','PaperPosition', [0 0 20 15])
set(gcf, 'InvertHardcopy', 'off')
set(gca,'color','black','FontName','Arial','FontSize',15,'LineWidth',2,'xcolor','white','ycolor','white','box','off')
% oa = oaxes
% 
% set(oa,'origin',[0 0 0],'ylabel',{'Attentional Modulation (dB power)' ''},'YLabelVerticalAlignment',{'middle'  'auto'})


axis(axes)
ylabel('Attentional Modulation (dB power)','FontName','Arial','FontSize',20)
xlabel('Time (s)','FontName','Arial','FontSize',20)
line([0 length(xData)],[0 0],'linewidth',1,'linestyle','--','color','w')
if lgnd
    legend(handles,lgnd,'FontName','Arial','FontSize',15,'textcolor','white','position',[0.81 0.88 0.0688 0.0547])
end
% line([xData(64) xData(181)], [0.12 0.12],'color',[.3 .6 1])
% line([xData(64) xData(64)], [0.12 0.11],'color',[.3 .6 1])
% line([xData(181) xData(181)], [0.12 0.11],'color',[.3 .6 1])
% text(xData(120),.125,'*','color',[.3 .6 1],'fontsize',15)

saveas(gcf,saveName,'fig')
saveas(gcf,saveName,'png')

end



