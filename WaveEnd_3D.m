function out = WaveEnd_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, timestore,start_indx, end_indx, time, saveon, Correlation)
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

    perc= 0.1;
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

%sometimes there are really weird pieces of the data that I chop out
%implicitly in the analysis but not here...so when we normalize, do so
%using smoothed data
casmooth = smoothdata(calcium); 
    cashort = [];
    calcium_demeaned = (calcium-min(casmooth))./(max(casmooth)-min(casmooth)); 

          
   for i = 1:length(end_indx)

    % run phase analysis for each oscillation
    cashort =  [calcium(start_indx(i):end_indx(i),:)];
    casmooth = smoothdata(cashort); 


           %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = (cashort - min(casmooth))./(max(casmooth)-min(casmooth));
 
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
       % MeanIslet= trimmean(calciumT',10);       % reference signal. 
       MeanIslet = mean(calciumT');
        %remove top and bottom 5% of data points just in case they are
        %super weird
        st = size(calciumT,1);
        for j=1:numcells % itterative index for cells
           [c(:,j)]=xcov(calciumT(:,j),MeanIslet,'none');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
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
       timehalf = find(calciumT(:,j)<0.5);
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
   [loc_mean, loc_spread, Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(length(indexes), ...
       Locations, newmaxCLvec(indexes,:), newmaxCLvec(indexes,:), cells_sorted(indexes,1:ceil(numcells*perc))', 1);



   % if saveon
   % csvwrite([savename '/' filename '/Wave_' 'Degree_KL.csv'], Degree_KL)
   % csvwrite([savename '/' filename '/Wave_' 'COG_KL.csv'], COG_KL)
   % end
  %% Quantify distance of all cells from center and radially: 
  IsletCenter = mean(Locations);
  dist_from_center = sqrt((Locations(:,1) - IsletCenter(1)).^2 +(Locations(:,2) - IsletCenter(2)).^2 + (Locations(:,3) - IsletCenter(3)).^2);
  dist_from_center = dist_from_center./max(dist_from_center); %normalize!




  %Check average duty cycle
    calcium_mean = rescale(mean(calcium(cuttime(i):cuttime(i+1),:), 2), 0, 1); 
    duty_cycle = length(find(calcium_mean > 0.5))./length(calcium_mean);%percent of time the oscillation is > 0.5
    if i == 1
        out.duty_cycle = duty_cycle;
    else
        out.duty_cycle_pka = duty_cycle;
    end

  %Check Frequency
  numosc = length(start_indx)./2;
    if i == 1
        freq = numosc./length(start_indx(1):cuttime(2)); %go until cuttime because sometimes there is a weird blimp after ading the solution
        out.freq = freq;
        out.loc_spread = loc_spread;
    else
        freq = numosc./length(start_indx(numosc+1):cuttime(3));
        out.freq_pka = freq;
        out.loc_spread_pka = loc_spread;
    end


  % NEED TO QUANTIFY RADIAL: 


    out.COG_KL_wave(i).trial = COG_KL;
    out.Degree_KL_wave(i).trial = Degree_KL;
        top10 = round(numcells.*perc);
   if i == 1
        out.MeanLocation = loc_mean;
        out.STLocation = loc_spread;
        out.dist_from_center = dist_from_center;
        out.top10 = cells_sorted(find(start_indx < cuttime(2)),1:top10);
        out.bottom10 = cells_sorted(find(start_indx < cuttime(2)), end-top10:end);
        for l = 1:length(start_indx)./2
            BottomLoc = [Locations(out.bottom10(l,:),1), Locations(out.bottom10(l,:),3), Locations(out.bottom10(l,:),3)];
            BottomCOG(l,:) = mean(BottomLoc);

            TopLoc = [Locations(out.top10(l,:),1), Locations(out.top10(l,:),3), Locations(out.top10(l,:),3)];
            TopCOG(l,:) = mean(TopLoc);
            out.Dist_traveled(l) = norm(TopCOG(l,:) - BottomCOG(l,:));
            out.phasespread_10(l) = (mean(newmaxCLvec(l,out.bottom10(l,:)))-mean(newmaxCLvec(l,out.top10(l,:))));
            out.Velocity(l) = out.Dist_traveled(l)./out.phasespread_10(l);
        end
        for k = 1:length(find(start_indx < cuttime(2)))
            for j = 1:length(find(start_indx < cuttime(2)))
                if k < j
                    intravar(k,j) = length(intersect(out.top10(k,:), out.top10(j,:)))/length(out.top10(k,:));
                    intravar_bottom(k,j) = length(intersect(out.bottom10(k,:), out.bottom10(j,:)))/length(out.bottom10(k,:));
                else
                    intravar(k,j) = NaN;
                    intravar_bottom(k,j) = NaN;
                end
            end
        end
        out.intravarall = intravar;
        out.intravar = mean2(intravar(~isnan(intravar)));
        out.intravar_bottom = mean2(intravar_bottom(~isnan(intravar_bottom)));
       for k = 1:size(out.top10,1)
            out.distcenter_top10(:,k) = dist_from_center(out.top10(k,:)); 
        end
        out.phase = newmaxCLvec(1:length(start_indx)./2,:);
        out.phasespread = range(newmaxCLvec(1:length(start_indx)./2,:)');
   else
       out.STLocation_pka = loc_spread;
       out.meanLocation_pka = loc_mean;        
       out.top10_pka = cells_sorted(find(start_indx > cuttime(2)),1:top10);
       out.bottom10_pka = cells_sorted(find(start_indx > cuttime(2)),end-top10:end);
    
        out.phase_pka = newmaxCLvec(length(start_indx)./2+1:end,:);
        out.dist_from_center_pka = dist_from_center;
        for k = 1:size(out.top10_pka,1)
            out.distcenter_top10_pka(:,k) = dist_from_center(out.top10_pka(k,:)) ;
            out.distcenter_bottom10_pka(:,k) = dist_from_center(out.bottom10_pka(k,:)) ;
        end
         for l = 1:length(start_indx)./2
            BottomLoc = [Locations(out.bottom10_pka(l,:),1), Locations(out.bottom10_pka(l,:),3), Locations(out.bottom10_pka(l,:),3)];
            BottomCOG(l,:) = mean(BottomLoc);

            TopLoc = [Locations(out.top10_pka(l,:),1), Locations(out.top10_pka(l,:),3), Locations(out.top10_pka(l,:),3)];
            TopCOG(l,:) = mean(TopLoc);
            out.Dist_traveled_pka(l) = norm(TopCOG(l,:) - BottomCOG(l,:));
            out.phasespread_10_pka(l) = (mean(newmaxCLvec(l+length(start_indx)./2,out.bottom10_pka(l,:)))-mean(newmaxCLvec(l+length(start_indx)./2,out.top10_pka(l,:))));
            out.Velocity_pka(l) = out.Dist_traveled_pka(l)./out.phasespread_10_pka(l);
        end

        
       for k = 1:length(find(start_indx > cuttime(2)))
            for j = 1:length(find(start_indx > cuttime(2)))
                if k < j
                    intravar_pka(k,j) = length(intersect(out.top10_pka(k,:), out.top10_pka(j,:)))/length(out.top10_pka(k,:));
                    intravar_bottom_pka(k,j) = length(intersect(out.bottom10_pka(k,:), out.bottom10_pka(j,:)))/length(out.bottom10_pka(k,:));

                else
                    intravar_pka(k,j) = NaN;
                    intravar_bottom_pka(k,j) = NaN;

                end
            end
       end
       out.intravar_pka_all = intravar_pka;
       out.intravar_pka = mean2(intravar_pka(~isnan(intravar_pka)));
       out.phasespread_pka = range(newmaxCLvec(length(start_indx)./2+1:end,:)');
       out.intravar_bottom_pka = mean2(intravar_bottom_pka(~isnan(intravar_bottom_pka)));
       
       for k = 1:length(find(start_indx > cuttime(2)))
            for j = 1:length(find(start_indx < cuttime(2)))
                intervar(k,j) = length(intersect(out.top10(j,:), out.top10_pka(k,:)))/length(out.top10_pka(k,:));
                intervar_bottom(k,j) = length(intersect(out.bottom10(j,:), out.bottom10_pka(k,:)))/length(out.bottom10_pka(k,:));

            end
       end
        out.intervar_all = intervar;
        out.intervar = mean2(intervar);
        out.intervar_bottom = mean2(intervar_bottom);
    end
    
  end

    [out.rot, out.rot_pka, out.rot_intra] = phase_poles(length(indexes), out.top10, out.bottom10, out.top10_pka, out.bottom10_pka, Locations);
     saveas(gcf, [savename '/' filename 'WaveAxis_repolarize.fig'])

       %spread of top and bottom: 
   
    %% find the average phase for all oscillations
    allphase = mean(finalphase);
        out.change = mean(abs(diff([mean(out.phase); mean(out.phase_pka)]./max(max([out.phase, out.phase_pka])))));
    


    [~,pre_top10] = sort(sum(out.phase));
    out.top10av = pre_top10(1:top10);
    [~,pka_top10] = sort(sum(out.phase_pka));
    out.top10pka = pka_top10(1:top10);
   out.top10percent_avg = length(intersect(out.top10av, out.top10pka))./top10;
     %check!!
     if find(size(out.phase) ~= size(out.phase_pka))
         error('The number of oscillations in pre and post are not the same')
     end


     if 0 %Turn off figures
     figure,
     for i = 1:wavenum
     nexttile
     plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(i,1:10)), 'linewidth',1, 'color', 'blue');
     hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(i,end-10:end)), 'linewidth',1, 'color', 'red');
    xline(time(round(start_indx(i)))), xline(time(round(end_indx(i))));
    xlim([time(round(start_indx(i))), time(round(end_indx(i)))]);
    title(['Oscillation Number ' num2str(i)])
    axg.data(i) = gca;
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})
    set(gca, 'color','none')
    xlabel('Time (s)')
    ylabel('Normalized Ca^{2+} Fluoresence')
    set(gca, 'box','off')

     end
    saveas(gcf, [savename '/' filename 'AllOscillations_repolarize.fig'])


% plot locations
    figure, 
   tt = tiledlayout(2, ceil(size(cells_sorted,1)./2));

    tt.TileSpacing = 'compact';
    tt.Padding = 'compact';
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

        saveas(gcf, [savename '/' filename '_repolarizeLocation.png'])
        saveas(gcf, [savename '/' filename '_repolarizeLocation.fig'])


    %% Calculate the number of high phase cells retained:
    
    %get top 10% of high phase cells
    top10 = round(numcells*.1);
    %first plot 'trajectory of top 5% of cells')
    topcells = cells_sorted(1,1:top10);
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
    retained(i) = length(intersect(cells_sorted(i,1:top10), topcells))./top10;
    end
    
    figure, bar(retained)
    ylabel('Percent of cells still in top 10%')
    xlabel('Oscillation')
        %saveas(gcf, [fileloc 'Highphasebar.png'])

   

  
        %saveAllFigsToPPT([savename '/' filename '/WaveInitiators_repolarize'])
     end
     close all
     

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