function [calcium, bo] = detrend_photobleaching(calcium)
%% for photobleached calcium oscillations
% Jennifer Briggs 11/2022

%% THIS PROGRAM FINDS A LINEAR REGRESSION OVER PEAKS OF CALCIUM OSCILLATIONS.
% assuming these peaks should be at the same height during second phase
% calcium oscillation, we can assume the linear decay in fluoresence is due
% to photobleaching and then remove. 

%check calcium direction: 
if size(calcium,1)>size(calcium,2)
    calcium = calcium';
end
[peak, yy] = findpeaks(mean(calcium), 'MinPeakProminence', 1*range(mean(calcium))/6);
[bo] = polyfit(yy,peak, 1);
%remove trendline
trendline = polyval(bo, [0:length(calcium)-1]);

%% if trend line isn't close to peaks: 
order_to_remove = [1, length(peak)];
i = 1;
peakhold = peak; yyhold = yy;
while max(abs(trendline(yy)'-peak')) > (max(mean(calcium))-min(mean(calcium)))/3
    peak = peakhold; yy = yyhold; 
    peak(order_to_remove(i)) = []; yy(order_to_remove(i)) = []; %remove some peaks
    [bo] = polyfit(yy,peak, 1);
    %remove trendline
    trendline = polyval(bo, [0:length(calcium)-1]);
    i = i+1
end

   
calcium = (calcium-trendline)./(trendline-min(min(calcium)));

