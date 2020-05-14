
function condOrder = pseudoRand(design,numTrials,dontRepeat)

%% function for pseudo-randomising condition orders
% design: EITHER
    % number of conditions (single integer)
    % m x n design matrix produced by fullfact (e.g. fullfact([2,3]), with m>1
% numTrials: number of trials in block/experiment
% dontRepeat: input 1 if do not want back-to-back repetitions of conditions


if nargin<3
    dontRepeat = 0;
end

% determine if input is number of conditions or design matrix
if size(design,1) > 1
    factorial = 1;
    numConds = length(design);
else 
    factorial = 0;
    numConds = design;
end

% allocate trial conditions
previous = nan;
condOrder = [];
condTally = zeros(1,numConds);
optionList = 1:numConds;
for T = 1:numTrials
    lowestTallyValue = min(condTally);
    if dontRepeat
        conditionOptions = find(condTally==lowestTallyValue & optionList~=previous);
    else
        conditionOptions = find(condTally==lowestTallyValue);
    end
    selectedID = randperm(length(conditionOptions),1);
    current = conditionOptions(selectedID);
    condOrder(T) = current;
    condTally(current) = condTally(current) + 1;
    previous = current;
end 


% if design matrix input assign conditions to each factor
if factorial
    for c = 1:numTrials
        for fact = 1:size(design,2)
            factors(fact,c) = design(condOrder(c),fact);
        end
    end
    condOrder = factors;
end




%% OLD VERSION

% if nargin<3
%     silent = 0;
% end
% 
% % determine if input is number of conditions or design matrix
% if size(design,1) > 1
%     factorial = 1;
%     numConds = length(design);
% else 
%     factorial = 0;
%     numConds = design;
% end
%     
% 
% % add equal number of trials for each condition
% numReps = floor(numTrials/numConds); % min number each condition
% condOrder=[];
% for c = 1:numConds
%     condOrder = [condOrder repmat(c,1,numReps)];
% end
% 
% % if unequal number of conditions, add random conditions
% addTrials = numTrials-length(condOrder);
% if addTrials
%     if ~silent
%         warning('Unequal number of conditions in design matrix')
%     end
%     addOrder = randperm(numConds);
%     condOrder = [condOrder addOrder(1:addTrials)];
% end
% 
% % randomise order
% order = randperm(numTrials);
% condOrder = condOrder(order);
% 
% % if design matrix input assign conditions to each factor
% if factorial
%     for c = 1:numTrials
%         for fact = 1:size(design,2)
%             factors(fact,c) = design(condOrder(c),fact);
%         end
%     end
%     condOrder = factors;
% end
% 


