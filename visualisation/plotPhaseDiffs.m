% function to plot phase angle differences 
%
% accepts vector of phase differences (e.g. between two electrodes at each
% timepoint)

function absAvg = plotPhaseDiffs(phaseDiffs)

avg = mean(exp(1i*phaseDiffs));
absAvg = abs(mean(exp(1i*phaseDiffs)));
THETA(2,:)=phaseDiffs'; % angles
RHO = repmat([0;1],1,length(THETA)); % set all lengths to 1
figure('position',[700 300 550 500])
polar(THETA,RHO)
hold on
l = polar([0;angle(avg)],[0; abs(avg)],'k');
set(l,'linewidth',5)

end
