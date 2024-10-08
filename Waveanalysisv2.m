%% by Jennifer Briggs 02/28/2022
%This script is modified off of code from Vira and Jenn - calculates phase of cells. 
%Input: calcium wave form
%Output: time of phase difference in ms
%V2 by EJ, can pick any number of waves and plot the top 10% of the cells 

    
%There are a few ways you can look at this wave form - either look at the
%whole thing or look at a single oscillation. Here we define what area we
%are look. Generally, if you are looking for a wave initiator, you'll want
%to look at just the beginning of an oscillation
    title('Select the beginning of 5 oscillations')
    starttime =  ginput(5)  %here you put the time in seconds that you want to start the analysis
    starttime = starttime(:,1);
    hold on, xline(starttime)
    title('Select the peak of those oscillations')
    endtime =    ginput(5)  %put the end time here
    endtime = endtime(:,1);
    xline(starttime, 'label','This is where we start')
    xline(endtime, 'label','This is where we end')
    wavenum=length(starttime);

    % for plotting: 
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
%% Start phase analysis MATLAB works by indexing the datapoints. Therefore, we
    %must find what index the time that you want to look at is.

    for i = 1:wavenum
    start_indx_f = find(abs(time-starttime(i))<0.5); 
    start_indx(i) = start_indx_f(1);
    end_indx_f = find(abs(time-endtime(i))<0.5);
    end_indx(i) = end_indx_f(1);
    end
    
    %maybe we should smooth the data;
   % calcium = smoothdata(calcium, 1, 'movmedian',3);

%% Wave origin analysis
    
    
%in order to make the analysis more accurate, we linearly interpolate between the
%points to artificially increase resolution. We may need to play with this
%as the outputs are not getting very good differentiation between phases.
    calcium_demeaned = (calcium-min(calcium))./(max(calcium)-min(calcium)); 
    %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = [];

        
          
   for i = 1:length(end_indx)
    % run phase analysis for each oscillation
    cashort =  [calcium_demeaned(start_indx(i):end_indx(i),:)];,...
       % calcium_demeaned(start_indx(2):end_indx(2),:);...
       % calcium_demeaned(start_indx(3):end_indx(3),:);...
       % calcium_demeaned(start_indx(4):end_indx(4),:);...
       % calcium_demeaned(start_indx(5):end_indx(5),:)];
    
 
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
       [c(:,j)]=xcov(calciumT(round(st/5):round(4*st/5),j),MeanIslet,'none');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
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
   
   
   [r_mean,Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(length(start_indx), Locations, newmaxCLvec(i,:), newmaxCLvec(i,:), phasevecsort(1:10));
   csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'Degree_KL.csv'], Degree_KL)
   csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'COG_KL.csv'], COG_KL)

   
    
   %% analyze the trajectory of cells
   for i = 1:numcells %identify the ranking for each cell in each oscillation
       [r, c] = find(cells_sorted' == i);
       ranking(:,i)  = 1-r./numcells; %If 1, first to depolarize
   end
   
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);
   figure, fig1 = heatmap(xvalues, yvalues, ranking')

   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')

   ylabel('Cell Number')
   xlabel('Oscillation Number')
   
   saveas(gcf, [fileloc 'HeatMap.png'])
    
  
    
    %% find the average phase for all oscillations
    allphase = mean(finalphase);
    
     figure,
     for i = 1:wavenum
     nexttile
     plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(i,1:(wavenum-1))), 'linewidth',1, 'color', 'blue')
     hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(i,end-3:end)), 'linewidth',1, 'color', 'red')
    xline(time(start_indx(i))), xline(time(end_indx(i)))
    title(['Oscillation Number ' num2str(i)])
    axg.data(i) = gca
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})
     end
          linkaxes([axg.data(1) axg.data(2) axg.data(3) axg.data(4) axg.data(5)], 'xy')

        saveas(gcf, [fileloc 'OscillationsWPhase.png'])
        saveas(gcf, [fileloc 'OscillationsWPhase.fig'])


     
    %X = Locations(:,1);
    %Y = Locations(:,2);
    %Z = Locations(:,3);
    %figure, scatter3(X,Y,Z,100, newmaxCLvec/max(newmaxCLvec), 'filled')
    %colormap hot
    %h = colorbar
    %set(gca, 'visible', 'off')
    
    %% Calculate the number of high phase cells retained:
    
    %get top 10% of high phase cells
    top10 = round(numcells*.1);
    %first plot 'trajectory of top 5% of cells')
    topcells = cells_sorted(1,1:top10)
    bottomcells = cells_sorted(1,end-top10+1:end);
    
    figure, plot(ranking(:, topcells),'k', 'linewidth',4)
    hold on, plot(ranking(:, bottomcells), 'b:','linewidth',3)
    
    %make legend
    ax =gca;
    axc = ax.Children;
    l1 = axc(end);
    l2 = axc(1);
    
    legend([l1, l2], '10% Highest Phase','10% Lowest Phase') 
    
    
    ylabel('Phase (1 = first to depolarize)')
    xlabel('Oscillation')
    saveas(gcf, [fileloc 'Highphasetraj.png'])
    
    %find percent of cells still within the top 5
    for i = 1:wavenum
    retained(i) = length(intersect(cells_sorted(i,1:top10), topcells))./top10
    end
    
    figure, bar(retained)
    ylabel('Percent of cells still in top 10%')
    xlabel('Oscillation')
        saveas(gcf, [fileloc 'Highphasebar.png'])

   

    
    
    
  
    
    
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

    for i = 1:5:length(calcium);

       nexttile(1)
       scatter3(X,Y,Z,100, calcium_demeaned(i,:), 'filled')
       %annotation('textarrow', mx,my,'String', 'Wave Initiators')
       text(mx+50,my-50,mz, 'Wave Initiators', 'FontSize',20)
       colorbar('south','axislocation', 'in')
       caxis manual
       caxis([0, 1])
       set(gca, 'visible', 'off')

       nexttile(2)
        ax = gca
        axc = ax.Children
        axc(1)
        axc(1).Value = i
        
       drawnow
%        if i == 1
%             keyboard
%         end
 
    end