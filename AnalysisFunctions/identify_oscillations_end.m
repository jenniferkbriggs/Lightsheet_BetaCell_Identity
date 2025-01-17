function [start_indx, end_indx] = identify_oscillations_end(start_indx, calcium, time)
%% Identify Oscillations manually or automatically
% Jennifer Briggs 11/2022

%use start index from end index in depolarization analysis 
if size(calcium,1)>size(calcium,2)
    calcium = calcium';
end

    wavenum = length(start_indx);
meancal = mean(calcium);
meandiff = diff(meancal);

figure, plot(meancal)

for i = 1:wavenum
    %select area around peak:
    start_indx_f = start_indx(i); 
    
    %find time of activity
    if start_indx_f(1)-200 < 1
        largearea = meancal(1:start_indx_f(1)+400);
    elseif start_indx_f(end)+450 > length(meancal)    
        largearea = meancal(start_indx_f(1)-200:end);
    else   
        largearea = meancal(start_indx_f(1)-200:start_indx_f(1)+400);
    end
    timeon = find((largearea - min(largearea))/(range(largearea)) > (0.3));
    timeon2 = find(diff(timeon) > 2);

    if ~isempty(timeon2)
        for j = 1:length(timeon2)
            if timeon2(j) < length(timeon)/2
                timeon = timeon(timeon2(j)+1:end);
            else
                try
                timeon = timeon(1:timeon2(j));
                catch
                    timeon = timeon(1:end-1);
                end
            end
        end
    end
    period = timeon(end)-timeon(1);
    
    st = start_indx_f(1)-round(period/4);
    if st < 1
        st = 1;
    end
    ed = start_indx_f(1) + round(5*period/4);
    
    meandiff2 = meandiff;
    meandiff2([1:st]) = 0;
    meandiff2([ed:end]) = 0;
    
   [pks,locs] = findpeaks(-(meancal(st:ed)-min(meancal(st:ed)))/range(meancal(st:ed)), 'MinPeakDist',10, 'MinPeakProminence',.2);


   if length(locs)<1 %then the window doesn't contain depolarization of the next oscillation
       end_indx(i) = ed;

   else
    end_indx(i) = locs(end)+st-period/10; %remember the location is indexed by a smaller window. Subtract by a bit because we don't want to be right up to the next depolarization
   end

    try
    if time(end_indx(i)) > starttime(i+1)
        disp('Error, oscillation is too thin to count')
        keyboard
    end
    end

    hold on, xline(start_indx(i), 'label',['start: ' num2str(i)]), xline(end_indx(i), 'color','r', 'label', ['end:' num2str(i)])
end
end