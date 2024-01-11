function out = Wave_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, timestore,start_indx, end_indx, time, saveon, Correlation)
%% by Jennifer Briggs 02/28/2022
%This script is modified off of code from Vira and Jenn - calculates phase of cells. 
%Input: calcium wave form
% ----: Location: number of cells x 3 (x,y,z) location
% ----: Correlation 1 if phase is determined using cross correlation and 0
%       if phase is determine using threshold
%Output: time of phase difference in ms
%V2 by EJ, can pick any number of waves and plot the top 10% of the cells 

    
%There are a few ways you can look at this wave form - either look at the
%whole thing or look at a single oscillation. Here we define what area we
%are look. Generally, if you are looking for a wave initiator, you'll want
%to look at just the beginning of an oscillation

%options:
if ~exist('saveon') %writing output to csv. If no options then write
    saveon = 1
end
    perc= 0.10
    wavenum=length(start_indx);
    numcells = size(calcium,2);

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
%% Wave origin analysis
    
    
%in order to make the analysis more accurate, we linearly interpolate between the
%points to artificially increase resolution. We may need to play with this
%as the outputs are not getting very good differentiation between phases.
if photobleaching
    [~,bo] = detrend_photobleaching(calcium(cuttime(1):cuttime(2),:)); % only remove photobleaching based on control!
    trendline = polyval(bo, [0:length(calcium)-1]');
    calstore = (calcium-trendline)./(trendline-min(min(calcium)));
else
    calstore = calcium;
end
    calcium_demeaned = (calcium-min(calcium))./(max(calcium)-min(calcium)); 
    %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = [];

        
          
   for i = 1:length(end_indx)
    % run phase analysis for each oscillation
    cashort =  [calcium_demeaned(start_indx(i):end_indx(i),:)];
       % calcium_demeaned(start_indx(2):end_indx(2),:);...
       % calcium_demeaned(start_indx(3):end_indx(3),:);...
       % calcium_demeaned(start_indx(4):end_indx(4),:);...
       % calcium_demeaned(start_indx(5):end_indx(5),:)];
    
 
    step = .005; %this gives how much to interpolate by
    xq = 1:step:size(cashort,1)+1;
    vq1 = interp1(1:size(cashort,1),cashort,xq); %may need to investigate a better way to do this
    %vq2 = spline(1:size(cashort,1),cashort',xq);
    
    timeunits = mean(diff(time), 'omitnan')*step*1000; %ms
    numcells=size(cashort,2);

    calciumT = (vq1);                           % new calcium time course
    [row,col] = find(isnan(calciumT)); %remove NaN's
    calciumT = calciumT(1:row(1)-1,:);
   % 
   % 
   %  % 2. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE
    clear vq1 cashort xq camax

    % 3. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET) WITH EACH INDIVIDUAL CELL
    tic

    %   HERE WE CALCULATE PHASE USING CROSS CORRELATIONS
    if Correlation
        MeanIslet= mean(calciumT,2);       % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
        st = size(calciumT,1);
        for j=1:numcells % itterative index for cells
           [c(:,j)]=xcov(calciumT(1:round(4*st/5),j),MeanIslet,'none');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
        end
        toc
    
        [maxCV, maxCL]=max(c);
        while length(unique(maxCL))<2;
            c(unique(maxCL), :) = [];
            [maxCV, maxCL]=max(c);
        end
        clear c
    else

   % HERE WE CALCULATE PHASE USING A THRESHOLD
   %make sure all cells are normalized: 
   calciumT = normalize(calciumT, "range");
   for j = 1:numcells
       timehalf = find(calciumT(:,j)>0.5);
       maxCL(j) = timehalf(1);
   end
  
    end
   % THE REST OF CODE IS USED FO

    newmaxCLvec_init = maxCL-mean(maxCL);
    newmaxCLvec(i,:) = newmaxCLvec_init/timeunits; %outputs phase difference in ms

   %  clear maxCL
    if isempty(nonzeros(newmaxCLvec_init))
        keyboard
        error('Something wrong with calcium signals')
    end
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
   
  % saveas(gcf, [fileloc 'HeatMap.png'])
  for i = 1:numtrial
    indexes = find(start_indx < cuttime(i+1) & start_indx > cuttime(i));
   [loc_mean, loc_spread, Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(length(indexes), Locations, newmaxCLvec(indexes,:), newmaxCLvec(indexes,:), cells_sorted(indexes,1:ceil(numcells*perc))', 1);
   % if saveon
   % csvwrite([savename '/' filename '/Wave_' 'Degree_KL.csv'], Degree_KL)
   % csvwrite([savename '/' filename '/Wave_' 'COG_KL.csv'], COG_KL)
   % end
  %% Quantify distance of all cells from center and radially: 
  IsletCenter = mean(Locations);
  dist_from_center = sqrt((Locations(:,1) - IsletCenter(1)).^2 +(Locations(:,2) - IsletCenter(2)).^2 + (Locations(:,3) - IsletCenter(3)).^2);

  % NEED TO QUANTIFY RADIAL: 


    
    out.COG_KL_wave(i).trial = COG_KL;
    out.Degree_KL_wave(i).trial = Degree_KL;
        top10 = round(numcells*.1);
   if i == 1
        out.MeanLocation = loc_mean;
        out.STLocation = loc_spread;
        out.top10 = cells_sorted(find(start_indx < cuttime(2)),1:top10);
        out.phase = newmaxCLvec;
        out.phaserange = max(newmaxCLvec, [],2)-min(newmaxCLvec, [],2); %ms
   else
       out.STLocation_pka = loc_spread;
       out.meanLocation_pka = loc_mean;        
       out.top10_pka = cells_sorted(find(start_indx > cuttime(2)),1:top10);
        out.phase_pka = newmaxCLvec;
        out.phaserange_pka = max(newmaxCLvec, [],2)-min(newmaxCLvec, [],2); %ms

        for i = 1:length(find(start_indx < cuttime(2)))
            for j = 1:length(find(start_indx < cuttime(2)))
                if i < j
                    intravar(i,j) = length(intersect(out.top10(i,:), out.top10(j,:)))/length(out.top10(i,:));
                else
                    intravar(i,j) = NaN;
                end
            end
        end
        out.intravarall = intravar;
        out.intravar = mean2(intravar(~isnan(intravar)));
        
       for i = 1:length(find(start_indx > cuttime(2)))
            for j = 1:length(find(start_indx > cuttime(2)))
                if i < j
                    intravar_pka(i,j) = length(intersect(out.top10_pka(i,:), out.top10_pka(j,:)))/length(out.top10_pka(i,:));
                else
                    intravar_pka(i,j) = NaN;
                end
            end
       end
       out.intravar_pka_all = intravar_pka;
       out.intravar_pka = mean2(intravar_pka(~isnan(intravar_pka)));
       
       
       for i = 1:length(find(start_indx > cuttime(2)))
            for j = 1:length(find(start_indx < cuttime(2)))
                intervar(i,j) = length(intersect(out.top10(j,:), out.top10_pka(i,:)))/length(out.top10_pka(i,:));
            end
       end
        out.intervar_all = intervar;
        out.intervar = mean2(intervar);
    end
    
  end
    %% find the average phase for all oscillations
    allphase = mean(finalphase);
    
     figure,
     for i = 1:wavenum
     nexttile
     plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(i,1:10)), 'linewidth',1, 'color', 'blue')
     hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(i,end-10:end)), 'linewidth',1, 'color', 'red')
    xline(time(start_indx(i))), xline(time(end_indx(i)))
    xlim([time(start_indx(i)), time(end_indx(i))])
    title(['Oscillation Number ' num2str(i)])
    axg.data(i) = gca
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})
    set(gca, 'color','none')
    xlabel('Time (s)')
    ylabel('Normalized Ca^{2+} Fluoresence')
    set(gca, 'box','off')

     end
     saveas(gcf, [savename '/' filename '/Wave_' 'AllOscillations.fig'])


% plot locations
    figure, 
   tt = tiledlayout(2, ceil(size(cells_sorted,1)./2))

    tt.TileSpacing = 'compact'
    tt.Padding = 'compact'
        for i = 1:size(cells_sorted, 1)
            nexttile
            scatter3(Locations(:,1), Locations(:,2), Locations(:,3), 75, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor',[0.7, 0.7, 0.7] , 'MarkerFaceAlpha', 0.5)
            hold on
            scatter3(Locations(cells_sorted(i,1:top10),1), Locations(cells_sorted(i,1:top10),2), Locations(cells_sorted(i,1:top10),3), 100, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue' )
            scatter3(Locations(cells_sorted(i,end-top10:end),1), Locations(cells_sorted(i,end-top10:end),2), Locations(cells_sorted(i,end-top10:end),3), 100, 'MarkerFaceColor', 'red', 'MarkerEdgeColor','red' )
            legend('Normal Cell','High Phase','Low Phase','Location', 'Northeast')
            set(gca, 'color','none')
    title(['Oscillation Number ' num2str(i)])

    xlabel('X (\mum)')
    ylabel('Y (\mum)')
    zlabel('Z (\mum)')
        end

        saveas(gcf, [savename '/' filename '/Wave_' 'HighPhaseLocation.png'])
        saveas(gcf, [savename '/' filename '/Wave_' 'HighPhaseLocation.fig'])


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
    %saveas(gcf, [fileloc 'Highphasetraj.png'])
    
    %find percent of cells still within the top 5
    for i = 1:wavenum
    retained(i) = length(intersect(cells_sorted(i,1:top10), topcells))./top10
    end
    
    figure, bar(retained)
    ylabel('Percent of cells still in top 10%')
    xlabel('Oscillation')
        %saveas(gcf, [fileloc 'Highphasebar.png'])

   

  
        saveAllFigsToPPT([savename '/' filename '/WaveInitiators'])

end

function out = watchoscillations(calcium_demeaned, cells_sorted, Locations)
     %% If you want to watch the calcium oscillations:
     addpath('~/Documents/GitHub/UniversalCode/gif/')
     X = Locations(:,1);
     Y = Locations(:,2);
     Z = Locations(:,3);
    figure, nexttile, 
    
    calcium_av= movmean(calcium_demeaned, 5);
    colormap 
 
    nexttile(2)
    plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
    hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(1:4)), 'linewidth',1, 'color', 'blue')
    hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(end-3:end)), 'linewidth',1, 'color', 'red')
    xline(1, 'linewidth', 4)
    set(gca, 'box','off')
    set(gcf, 'color','white')
    xlabel('Time (s)')
    ylabel('Calcium Fluoresence')
    


    xwaveinit = X(cells_sorted(1:5));
    ywaveinit = Y(cells_sorted(1:5));
    zwaveinit = Z(cells_sorted(1:5));
    
    my = mean(ywaveinit);
    mx = mean(xwaveinit);
    mz = mean(zwaveinit);

    gif('~/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/Results/CalciumWave.gif')
    for i = 1:3:length(calcium);

       nexttile(1)
       scatter3(X,Y,Z,100, calcium_demeaned(i,:), 'filled')
       %annotation('textarrow', mx,my,'String', 'Wave Initiators')
       text(mx+10,my-10,mz, 'Wave Initiators', 'FontSize',20)
       c = colorbar('south','axislocation', 'in')
       ylabel(c, 'Phase')
       caxis manual
       caxis([0, 0.8])
       ax = gca
       ax.XTick = [];
       ax.YTick = [];
       ax.ZTick = [];


       nexttile(2)

     
        ax = gca;
        axc = ax.Children;
        axc(1)
        axc(1).Value = time(i)
        
       drawnow
%        if i == 1
%             keyboard
%         end
    gif
 
    end
end