% function yfilt = gaussFilter1D(y,sigma,size)
% 
% sigma = 3;
% size = 15;
% x = linspace(-size / 2, size / 2, size);
% gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
% gaussFilter = gaussFilter / sum (gaussFilter); % normalize
% 
% 
% % y = rand(500,1);
% 
% % yfilt = filter (gaussFilter,1, y);
% yfilt = conv (y, gaussFilter, 'same');

sigma = 400
width = length(tone)/2;
support = (-width:width);
gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
gaussFilter = gaussFilter/ sum(gaussFilter);

close all
% yfilt = conv (y, gaussFilter, 'same');
% yfilt = gaussFilter.*tone
gaussFilter = gaussFilter(1:length(tone))
% plot(gaussFilter)

toneFilt = tone.*gaussFilter
plot(toneFilt)