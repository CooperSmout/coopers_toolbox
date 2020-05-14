
function [data, classList, numReps, savedTrialIDX, groupOrder] = groupTrialsByClass(data,classList,numTrials2avg)

    % identify categories
    categories = unique(classList);
    categories = categories(~isnan(categories)); % remove nans
    numCategories = length(categories);
    
    % tally how many trials in each category
    categoryTally = zeros(1,numCategories);
    for C = 1:numCategories
        categoryTally(C) = sum(classList==categories(C)); % tally orientations 
    end
    
    % work out how many groups there will be for each category 
    numReps = floor(min(categoryTally)/numTrials2avg);
    
    % allocate groups and average trials
    groupOrder = zeros(numCategories,numReps,numTrials2avg);
    savedTrials = false(1,length(classList));
    for C = 1:numCategories
        catX = find(classList==categories(C));
        keepX = catX(randperm(length(catX),numReps*numTrials2avg)); % trim unused trials
        groupOrder(C,:,:) = sort(reshape(keepX,numReps,numTrials2avg)); 
        for rep = 1:numReps
            saveThisTrial = groupOrder(C,rep,1);
            savedTrials(saveThisTrial) = 1;
            data(:,:,saveThisTrial) = mean(data(:,:,groupOrder(C,rep,:)),3); % collapse trials
        end
    end
    
    % remove collapsed trials
    data(:,:,~savedTrials) = [];
    classList(~savedTrials) = [];
    
    % record which trials actually kept
    savedTrialIDX = find(savedTrials);
    
    
end
    
    
%     % OLD
%     % group (average) trials
%     tally = 1;
%     groupedData = nan(size(data,1),size(data,2),numReps*numCategories);
%     groupedCategories = nan(1,numReps*numCategories);
%     for C = 1:numCategories
%         currentID = find(classList==cats(C));
%         currentID = currentID(randperm(length(currentID),numReps*numTrials2avg)); % trim unused trials and randomise order
%         for rep = 1 : numTrials2avg : numReps*numTrials2avg
%             groupedData(:,:,tally) = mean(data(:,:,currentID(rep:rep+numTrials2avg-1)),3);
%             groupedCategories(tally)     = C;
%             tally = tally + 1;
%         end
%     end
%     
%     