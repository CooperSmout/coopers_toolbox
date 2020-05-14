

function plotConnectivityByDistance(dat1,dat2,elecs)

% variables
if nargin<3
    elecs = 1:64;
end



% plot connectivity by distance
load biosemi64.mat
dat=[];
for E1 = 1:64
    for E2 = 1:E1-1;
        
        if any(elecs==E1) && any(elecs==E2)
        
            dat(end+1,1)    = biosemi64.chanDistances(E1,E2);
            dat(end,2)      = mean(dat1(E1,E2,:),3);
            dat(end,3)      = mean(dat2(E1,E2,:),3);
        end
    end
end
figure; hold on
scatter(dat(:,1),dat(:,2),'r.')
scatter(dat(:,1),dat(:,3),'k.')
ylim ([0 1])
legend({inputname(1),inputname(2)})
title('Phase Coherence by Electrode Separation')
xlabel('electrode separation')
ylabel('phase coherence')


% plot mean line
l1=line([min(dat(:,1)) max(dat(:,1))],[mean(dat(:,2)) mean(dat(:,2))]);
set(l1,'color',[1 0 0],'linestyle','--')
l2=line([min(dat(:,1)) max(dat(:,1))],[mean(dat(:,3)) mean(dat(:,3))]);
set(l2,'color',[0 0 0],'linestyle','--')


