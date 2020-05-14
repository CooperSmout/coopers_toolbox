%% function for pseudo-randomising condition orders

% design: EITHER
    % number of conditions (single integer)
    % design matrix produced by fullfact (e.g. fullfact([2,3])
% numTrials: number of trials in block/experiment

function condOrder = pseudoRandReplace(numConds,numTrials)


for level = 1:2
    condOrder(level,:) = pseudoRand(numConds,numTrials);
    if level>1
        while any(condOrder(level,:)==condOrder(level-1,:))
            condOrder(level,:) = pseudoRand(numConds,numTrials);
        end
    end
end


