function [start_indx, end_indx] = identify_oscillations(calcium, time, auto)
%% Identify Oscillations manually or automatically
% Jennifer Briggs 11/2022
if size(calcium,1)>size(calcium,2)
    calcium = calcium';
end
if auto ==0
    title('Select approximate peak of interesting oscillations')
    starttime =  ginput()  %here you put the time in seconds that you want to start the analysis
    starttime = starttime(:,1);
else   
    cal = (mean(calcium')-min(mean(calcium')))./range(mean(calcium'));
    [~, starttime] = findpeaks(cal, time, 'MinPeakProminence', .05, 'MinPeakDist',60);
    starttime(starttime < 200) = [];
end

wavenum = length(starttime);

    for i = 1:wavenum
    xline(starttime(i), 'label',['wavenumber: ' num2str(i)])
    end
    
    title('Click continue or type dbcont if happy')
    
%keyboard
    
meancal = mean(calcium);
meandiff = diff(meancal);

figure, plot(meancal)

for i = 1:wavenum
    %select area around peak:
    start_indx_f = find(abs(time-starttime(i))<0.5); 
    
    %find time of activity
    if start_indx_f(1)-200 < 1
        largearea = meancal(1:start_indx_f(1)+200);
    elseif start_indx_f(end)+200 > length(meancal)    
        largearea = meancal(start_indx_f(1)-200:end);
    else   
        largearea = meancal(start_indx_f(1)-200:start_indx_f(1)+200);
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
    
    st = start_indx_f(1) - round(2*period/3);
    if st < 1
        st = 1;
    end
    ed = start_indx_f(1) + round(period/4);
    
    meandiff2 = meandiff;
    meandiff2([1:st]) = 0;
    meandiff2([ed:end]) = 0;
    
    stinx = find(meandiff2 == max(meandiff2(st:ed)));
    stinx = stinx(stinx > st & stinx < ed)-round(period/4);

    start_indx(i) = stinx(1);
    
  
    [pks,locs] = findpeaks((meancal(st:ed)-min(meancal(st:ed)))/range(meancal(st:ed)), 'MinPeakDist',10, 'MinPeakProminence',.2);
    if isempty(locs)
         [pks,locs] = max(meancal(st:ed));
    end
    [~,minpkloc] = min(abs(start_indx_f(1)-locs));
    edinx = max([st+100, st+locs(minpkloc)-period/8]); %%There is a minimum length that must be evaluated over
       
    edinx = edinx(edinx > st & edinx < ed+100);
    end_indx(i) = round(edinx(1));
    try
    if time(end_indx(i)) > starttime(i+1)
        disp('Error, oscillation is too thin to count')
        keyboard
    end
    end

    hold on, xline(start_indx(i), 'label',['start: ' num2str(i)]), xline(end_indx(i), 'label', ['end:' num2str(i)])
end
end