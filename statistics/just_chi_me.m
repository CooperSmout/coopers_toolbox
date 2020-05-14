

% calculate frequencies of synchronous judgements and plot
[famTable,CHI2,P,labels] = crosstab(famResults(:,1),famResults(:,2));
famHist = famTable(:,2)./sum(famTable,2);
bar(famHist)
title(['Practice ' opt.participant '   p=' num2str(P)])


% calculate frequencies of synchronous judgements and plot
[mcsTable, CHI2, P, labels] = crosstab(mcsResults(:,1),mcsResults(:,2));
mcsHist = mcsTable(:,2)./sum(mcsTable,2);
bar(mcsHist)
title(['MCS ' opt.participant '   p=' num2str(P)])



% [~, famP]=chiSquareTest(mcsPractice,.1)
% mcsPracticeHist = mcsPractice(2,:)./sum(mcsPractice)
% bar(mcsPracticeHist)
% ylim([0 1])
% title(['FAM participant ' opt.participant ', p=' num2str(famP)])
% saveas(gcf,[ opt.participant ' fam'],'jpg')
% 
% bar(mcsTestHist)
% ylim([0 1])
% title(['MCS participant ' opt.participant ', chi-square p=' num2str(P)])
% saveas(gcf,[ opt.participant ' mcs'],'jpg')
