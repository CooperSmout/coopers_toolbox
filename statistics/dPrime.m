%% Calculate percent correct and d' for each coherency value

function PC = dPrime(data)

PC = unique(data(:,1));

for i=1:length(PC)

    id = (data(:,1) == PC(i)) & (data(:,2) ~= NaN);
    PC (i,2) = sum(data(id,2) == 1);
    PC (i,3) = sum(id);
    PC (i,4) = PC(i,2)/PC(i,3); % divide no. correct responses by total number of responses
    PC (i,5) = sqrt(2)*norminv(PC(i,4)); % estimate d' for each value of coherency
    
%     PC (i,4) = PC (i,2)*100; % multiply into percentage for display only
    
end

% disp(['number trials so far = ' num2str(length(data))])


end