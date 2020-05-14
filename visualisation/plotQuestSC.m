function plotQuestSC(RF,PF)

error('deprecated, use PAL_plotStaircase')

questParams = [RF.mean RF.beta RF.gamma RF.lambda];

% plot staircase 
figure('Position', [400 300 600 500]);
hold on
plot(1:100,PF(questParams,1:100),'b')

t = 1:length(RF.x);
figure('name','Running Fit Adaptive Procedure','Position', [500 400 600 500]);
plot(t,RF.x,'k');
hold on;
plot(t(RF.response == 1),RF.x(RF.response == 1),'ko', ...
    'MarkerFaceColor','k');
plot(t(RF.response == 0),RF.x(RF.response == 0),'ko', ...
    'MarkerFaceColor','w');
set(gca,'FontSize',16);
axis([0 length(RF.x)+1 min(RF.priorAlphaRange) max(RF.priorAlphaRange)])
line([1 length(RF.x)],[RF.mean RF.mean],'linewidth', 1,'linestyle', '--', 'color','b');
xlabel('Trial');
ylabel('Stimulus Intensity');