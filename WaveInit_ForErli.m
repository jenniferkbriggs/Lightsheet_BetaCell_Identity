%% by Jennifer Briggs 02/28/2022
%This script is modified off of code from Vira and Jenn - calculates phase of cells. 
%Input: calcium wave form
%Output: time of phase difference in ms

clear all
close all
clc

w = warning('off','all');

try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

%change the percent of cells you are looking at
   percent_to_analyze = 0.05
   
%decide whether to automatically or manually pick waves
    auto = 1; % automatically pick oscillations

%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)

fileloc = ('/Volumes/Briggs_10TB/Merrin/Ca_Courses/Singlecelltraces/EJ106/')
cd(fileloc)
files = dir('.*pos*')
for fi = 1:length(files)
filename = files(fi).name(3:9)
fullfile = [fileloc filename]
%Here we load the file
    calcium = readmatrix([fullfile '_Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    figure, plot(time,calcium) %plot calcium

    disp("Resize the figure and then click continue after you are happy with it")
    %keyboard
    Locations = readmatrix([fullfile ' pos_Detailed.csv']);
    Locations = Locations(:,1:3);
    
%There are a few ways you can look at this wave form - either look at the
%whole thing or look at a single oscillation. Here we define what area we
%are look. Generally, if you are looking for a wave initiator, you'll want
%to look at just the beginning of an oscillation

calcium = smoothdata(calcium, 1, 'movmean',20);
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
    title(filename)
%% Start phase analysis MATLAB works by indexing the datapoints. Therefore, we
    %must find what index the time that you want to look at is.
% 
meancal = mean(calcium');
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
    edinx = st+locs(minpkloc)-period/8;
    edinx = edinx(edinx > st & edinx < ed+100);
    end_indx(i) = round(edinx(1));
    
    hold on, xline(start_indx(i), 'label',['start: ' num2str(i)]), xline(end_indx(i), 'label', ['end:' num2str(i)])
    
    

end
    

%% Wave origin analysis
    
    
%in order to make the analysis more accurate, we linearly interpolate between the
%points to artificially increase resolution. We may need to play with this
%as the outputs are not getting very good differentiation between phases.
    calcium_demeaned = calcium;
    %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = [];

        
          
   for i = 1:length(end_indx)
    % run phase analysis for each oscillation
    cashort =  [calcium_demeaned(start_indx(i):end_indx(i),:)];
    
    %demean
    cashort = (cashort-min(cashort))./(max(cashort)-min(cashort));
 
    step = .01; %this gives how much to interpolate by
    xq = 1:step:size(cashort,1)+1;
    vq1 = interp1(1:size(cashort,1),cashort,xq); %may need to investigate a better way to do this
    %vq2 = spline(1:size(cashort,1),cashort',xq);
    
    timeunits = mean(diff(time))*step*1000; %ms
    numcells=size(cashort,2);

    calciumT = (vq1);                           % new calcium time course
    [row,col] = find(isnan(calciumT)); %remove NaN's
    calciumT = calciumT(1:row(1)-1,:);
    % 2. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE
    MeanIslet= mean(calciumT,2);       % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
    clear vq1 cashort xq camax

    % 3. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET) WITH EACH INDIVIDUAL CELL
    tic
    st = size(calciumT,1);
    for j=1:numcells % itterative index for cells
       [c(:,j)]=xcov(calciumT(round(st/5):round(4*st/5),j), MeanIslet, 'none');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    end
    toc

    [maxCV, maxCL]=max(c);
    clear c

   % 4. PLOTTING SIGNAL, XCOV, AND OUTPUTTING MAX XCOV AND CORRESPONDING TIME LAG
    newmaxCLvec_init = maxCL-mean(maxCL);
    newmaxCLvec(i,:) = newmaxCLvec_init/timeunits; %outputs phase difference in ms


    %this is where you get the final output. phasevecsort gives you the
    %sorted vector of the phase lag compared to the islet mean.
    %cells_sorted is what you are really interested in. This gives the cell
    %index in order from high phase (start earlier) to low
    %phase (start later)
    [phasevecsort_init, cells_sortedinit] = sort(newmaxCLvec(i,:)); 
    cells_sorted(i,:)= cells_sortedinit;
    phasevecsort(i,:) = phasevecsort_init;
    finalphase = newmaxCLvec;
    
   end
   
   %rank cells within oscillation
   for i = 1:numcells %identify the ranking for each cell in each oscillation
       [r, c] = find(cells_sorted' == i);
       ranking(:,i)  = 1-r./numcells; %If 1, first to depolarize
   end
   
   %% Analysis: 
   
   opts.figs = 0;
   % look at the location of high phase and low phase cells:
   for i = 1:wavenum
     [high_dist, low_dist, pos_new, highphasecenter(i,:), lowphasecenter(i,:), V, D]  = locmap(Locations, ranking(i,:), percent_to_analyze, opts); %can play with last number, this 
      dist_from_high(:,i) = high_dist;
      dist_from_low(:,i) = low_dist;
      loc_along_wave(:,i) = pos_new;
      eigenvec(:,i) = V(:,1);
      eigenval(i) = D(1,1);
   end
   
   
%visualize iset
figure, scatter3(Locations(:,1), Locations(:,2), Locations(:,3), 40, [0.8, 0.8, 0.8])
ax = gca
ax.GridLineStyle = 'none'
oscillationcolor = jet(wavenum);
hold on, s= scatter3(highphasecenter(:,1), highphasecenter(:,2), highphasecenter(:,3), 100,oscillationcolor, '*')
hold on, e=scatter3(lowphasecenter(:,1), lowphasecenter(:,2), lowphasecenter(:,3), 100, oscillationcolor, 'filled')
h = colorbar;
ylabel(h, 'Oscillation Number')
h.Ticks = [0:1/(wavenum-1):1]
h.TickLabels = sprintfc('%d',[1:wavenum])
legend([s(1) e(1)], 'High Phase Center', 'Low Phase Center', 'location','east')

%plotting things:
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);

   ydisplayvalues = sprintfc('%d',[1:(numcells)]);
  ydisplayvalues(1:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(2:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(3:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(4:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(5:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(6:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(7:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(8:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(9:10:length(ydisplayvalues)) = {''};

   

   %disp: 
   %sort based on average 
   [~, dist_from_high_sort] = sort(mean(dist_from_high,2));
   figure,  heatmap(xvalues, yvalues, dist_from_high(dist_from_high_sort,:))
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   ylabel('Ranking sorted by average value')
   xlabel('Oscillation Number')
    title('Cell distance (normalized by diameter) from center of highest phase')

    
   [~, dist_from_low_sort] = sort(mean(dist_from_low,2));
   figure,  heatmap(xvalues, yvalues, dist_from_low(dist_from_low_sort,:))
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   ylabel('Ranking sorted by average value')
   xlabel('Oscillation Number')
    title('Cell distance (normalized by diameter) from center of lowest phase')
   
   [~, loc_along_wave_sort] = sort(mean(loc_along_wave,2));
   figure,  heatmap(xvalues, yvalues, loc_along_wave(loc_along_wave_sort,:))
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   ylabel('Cell Number')
   ylabel('Ranking sorted by average value')
    title('Cell location along wave')
   
   
%% analyze the trajectory of cells

  
   [~, ranking_sort] = sort(mean(ranking,1));

   figure, fig1 = heatmap(xvalues, yvalues, ranking(:,ranking_sort)');
   ydisplayvalues = sprintfc('%d',[1:(numcells)]);
  ydisplayvalues(1:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(2:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(3:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(4:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(5:10:length(ydisplayvalues)) = {''};
      ydisplayvalues(6:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(7:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(8:10:length(ydisplayvalues)) = {''};
   ydisplayvalues(9:10:length(ydisplayvalues)) = {''};

   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')

   ylabel('Ranking sorted by average value')
   xlabel('Oscillation Number')
   title('Cell Phase')
   %saveas(gcf, [fileloc 'HeatMap.png'])
    
  
    
    %% find the average phase for all oscillations
    allphase = mean(finalphase);
    
   calcium_demeaned = (calcium_demeaned - min(calcium_demeaned))./(range(calcium_demeaned));

    
     figure,
     for i = 1:wavenum
     nexttile
     plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(i,1:length(ranking)*percent_to_analyze)), 'linewidth',1, 'color', 'blue');
    hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(i,end-percent_to_analyze*length(ranking):end)), 'linewidth',1, 'color', 'red');
    xline(time(start_indx(i))), xline(time(end_indx(i)));
    xlim([time(start_indx(i)), time(end_indx(i))])
    title(['Oscillation Number ' num2str(i)])
    axg.data(i) = gca;
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})
     end
        %  linkaxes([axg.data(1) axg.data(2) axg.data(3) axg.data(4) axg.data(5)], 'xy')

   %     saveas(gcf, [fileloc 'OscillationsWPhase.png'])
   %     saveas(gcf, [fileloc 'OscillationsWPhase.fig'])

     
        
    X = Locations(:,1);
    Y = Locations(:,2);
    Z = Locations(:,3);
    figure,
    for i = 1:wavenum
     nexttile,
     scatter3(X,Y,Z,100, newmaxCLvec(i,:)/max(newmaxCLvec(i,:)), 'filled')
    colormap hot
    h = colorbar;
    set(gca, 'visible', 'off')
    title(['Oscillation number ' num2str(i)])
    end
    %% Calculate the number of high phase cells retained:
    
    %get top 5% of high phase cells
    top5 = round(numcells*.1);
    %first plot 'trajectory of top 5% of cells')
    topcells = cells_sorted(1,1:top5);
    bottomcells = cells_sorted(1,end-top5+1:end);
    
    figure, plot(ranking(:, topcells),'k', 'linewidth',4)
    hold on, plot(ranking(:, bottomcells), 'b:','linewidth',3)
    
    %make legend
    ax =gca;
    axc = ax.Children;
    l1 = axc(end);
    l2 = axc(1);
    
    legend([l1, l2], [num2str(percent_to_analyze*100) '%  Highest Phase'], [num2str(percent_to_analyze*100) '%  Lowest Phase']) 
    
    
    ylabel('Phase (1 = first to depolarize)')
    xlabel('Oscillation')
    saveas(gcf, [fileloc 'Highphasetraj.png'])
    
    %find percent of cells still within the top 5
    for i = 1:wavenum
        for j = 1:wavenum
    retained(i,j) = length(intersect(cells_sorted(i,1:top5), cells_sorted(j,1:top5)))./top5;
        end
    end
    
    figure, 
    for i = 1:wavenum
   nexttile, bar(retained(i,:))
    ylabel(['Percent of cells still in top ' num2str(percent_to_analyze*100) '%'])
    xlabel('Oscillation')
    title(['Oscillation: ' num2str(i)])
    end
      %  saveas(gcf, [fileloc 'Highphasebar.png'])
      clearvars -except   percent_to_analyze auto fileloc files
    try
   saveAllFigsToPPT([fullfile 'PhaseAnalysis'])
  end
end
    
    
    %% If you want to watch the calcium oscillations:
    figure, nexttile, 
    
    calcium_av= movmean(calcium_demeaned, 5);
    colormap 
 
    nexttile(2)
    plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
    hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(1:4)), 'linewidth',1, 'color', 'blue')
    hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(end-3:end)), 'linewidth',1, 'color', 'red')
    xline(1, 'linewidth', 4)

    xwaveinit = X(cells_sorted(1:5));
    ywaveinit = Y(cells_sorted(1:5));
    zwaveinit = Z(cells_sorted(1:5));
    
    my = mean(ywaveinit);
    mx = mean(xwaveinit);
    mz = mean(zwaveinit);
    
    
% 
%     for i = 1:5:length(calcium);
% 
%        nexttile(1)
%        scatter3(X,Y,Z,100, calcium_demeaned(i,:), 'filled')
%        %annotation('textarrow', mx,my,'String', 'Wave Initiators')
%        text(mx+50,my-50,mz, 'Wave Initiators', 'FontSize',20)
%        colorbar('south','axislocation', 'in')
%        caxis manual
%        caxis([0, 1])
%        set(gca, 'visible', 'off')
% 
%        nexttile(2)
%         ax = gca
%         axc = ax.Children
%         axc(1)
%         axc(1).Value = i
%         
%        drawnow
% %        if i == 1
% %             keyboard
% %         end
%  
%     end