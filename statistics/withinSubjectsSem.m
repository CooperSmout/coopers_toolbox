
% calculates means and 95% confidence intervals for within subjects design
% (for plotting purposes only)
% see Cousineau (2005) and Morey (2008)

% input data should have participants in rows and within-subject conditions in columns
% accepts second input of 'CI' to output confidence interval instead of
% standard error of the mean


function [semOrCI, m] = withinSubjectsSem(data,varargin)

    
% calculate means
m = mean(data);


% normalisation method proposed by Cousineau (2005)
subjMeans = mean(data,2); % calculate mean for each participant
for cond = 1:size(data,2)
    meanCentredData(:,cond) = data(:,cond) - subjMeans; % subtract each participants mean from each condition score
end
newData = meanCentredData + mean(data(:)); % add grand mean to each score

cousineauSem = std(newData)/sqrt(size(newData,1));
cousineauSem = cousineauSem(1); % because all values same across conditions
cousineauCI = 1.96 * cousineauSem; % 95% confidence interval according to Cousineau (2005)


% correction to Cousineau calculation (Morey, 2008; Franz, 2012)
moreySem = cousineauSem * sqrt( size(data,2)/(size(data,2)-1) ) ; 
moreyCI = 1.96 * moreySem; 


% determine if sem or CI output
if nargin==1
    semOrCI = moreySem;
elseif strcmp(varargin{1},'CI')
    semOrCI = moreyCI;
end


end

