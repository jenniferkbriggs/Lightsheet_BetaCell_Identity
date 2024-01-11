function out = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, timestore,start_indx, end_indx, saveon)
%% Run Network Analysis
% Jennifer Briggs 06/02/2022
ct = 1;
%options:
if ~exist('saveon') %writing output to csv. If no options then write
    saveon = 1
end

if photobleaching
    [~,bo] = detrend_photobleaching(calcium(cuttime(1):cuttime(2),:)); % only remove photobleaching based on control!
    trendline = polyval(bo, [0:length(calcium)-1]');
    calstore = (calcium-trendline)./(trendline-min(min(calcium)));
else
    calstore = calcium;
end
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);
cellnum = min([size(calcium,1), size(calcium,2)]);
Opts.Max = cellnum/8;
Opts.Method = 'Degree'
Opts.avDeg = 22;
Threshold(i) = findoptRth(calcium, Opts);
end


%Threshold = mean(Threshold) %outstanding question - should I take the average of the thresholds??
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

Opts.direction_hold = 1;
Opts.printasSubplot = 0;
Opts.Subplotnum = 0;
Opts.figs = 1;
Opts.multiseed = 0;
Opts.multiseednum = 1;
fignum = 1

disp(Threshold(i))
[N, adj, kperc, histArrayPercShort, Rij] = links(calcium,Threshold(i),Opts,fignum); %%This is where the network is built
[sorted, cellsor]= sort(N);
%network analysis is performed
[L, EGlob, CClosed, ELocClosed, COpen, ELocOpen, nopath]  = graphProperties(adj);
G = graph(adj, 'upper');
Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);

%%------plot network------%%
try
    c = 0.7.*ones(size(Locations));
    c(Hubs,1) = ones(length(Hubs), 1).*136/255;
    c(Hubs,2) = ones(length(Hubs), 1).*8./255';
    c(Hubs,3) = ones(length(Hubs), 1).*8./255';
    figure, plot(graph(sparse(adj), 'upper'), 'Xdata',Locations(:,1), 'Ydata',Locations(:,2), ...
        'Zdata',Locations(:,3),'NodeColor',c, 'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3)
    %hidden graphs for legend
    hold on,
    ghubs = plot(graph(sparse([0,1;1,0])), 'NodeColor',[136/255,8./255,8./255'],...
        'Xdata',Locations(1:2,1),'Ydata',Locations(1:2,2),'Zdata',Locations(1:2,3),...
        'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3)
    hold on,
    gnonhubs = plot(graph(sparse([0,1;1,0])), 'NodeColor',[0.7,0.7,0.7],...
        'Xdata',Locations(1:2,1),'Ydata',Locations(1:2,2),'Zdata',Locations(1:2,3),...
        'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3)
    legend([ghubs, gnonhubs], {'Hubs','Non-hubs'},'Location', 'northeast')
    set(gca, 'color','none')
    set(gca, 'box','off')
    set(gca,'ytick',[]), set(gca,'xtick',[]), set(gca,'ztick',[])
    
    %turn off graphs for legend
    set(ghubs,'Visible','off'), set(gnonhubs,'Visible','off')
end
    top10 = length(N)*.1;
    
if i == 1
    for j = 1:5
     figure(j)
     ax = gca
     if j<6
        text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'Control','FontSize',30)
     else
        text(ax.XLim(end)-range(ax.XLim)/4, ax.ZLim(end)-range(ax.ZLim)/4,'Control','FontSize',30)
     end
    end
    Hubs_control = Hubs;
    Netinfo_control = [L, EGlob, CClosed,ELocOpen, nopath, mean2(Rij)];
    out.Hubs = Hubs;
    out.L = L;
    out.Threshold_all = Threshold;
    out.N = N;
    out.Corr = mean2(Rij);
    out.sorted_deg = flipud(cellsor);
else
    for j = 1:5
        figure(j+5)
        ax = gca
        if j <6
            text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'PKa Application','FontSize',30)
        else
            text(ax.XLim(end)-range(ax.XLim)/4, ax.ZLim(end)-range(ax.ZLim)/4,'PKa Application','FontSize',30)
        end
    end
    Hubs_pka = Hubs;
    Netinfo_pka = [L, EGlob, CClosed,ELocOpen, nopath, mean2(Rij)];
    out.Hubs_pka = Hubs;
    out.L_pka = L;
    out.N_pka = N;
    out.Corr_pka = mean2(Rij);
    out.sorted_deg_pka = flipud(cellsor);
    out.mainted_perc = length(intersect(out.sorted_deg(1:top10), out.sorted_deg_pka(1:top10)))/length(out.sorted_deg(1:top10))
end


end
%percent of Hubs that are consistent across conditions:



    
%% analysis for consistency in hubs across oscillations:
% for plotting: 
numcells = length(adj)
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
       
start_indx_hold = start_indx;

for i = 1:numtrial
    calcium = calstore;
    time = timestore(cuttime(i):cuttime(i+1),:);
      if size(calcium,1)<size(calcium,2)
         calcium = calcium';
      end

    %%
    calcium = (calcium-min(calcium))./range(calcium); %demean
    if i == 1 & numtrial == 2
        start_indx = start_indx_hold(start_indx_hold < cuttime(i+1) & start_indx_hold >= cuttime(i));
        next_start_indx = find(start_indx_hold < cuttime(i+1)); %defines the end of the cycle
        start_indx= [start_indx cuttime(2)];
    else
        start_indx = start_indx_hold(start_indx_hold < cuttime(i+1) & start_indx_hold >= cuttime(i));
        start_indx  = [start_indx cuttime(end)];
    end
    
    for j = 1:length(start_indx)-1 %the last start_indx is the end of the trial
        if start_indx(j) < 1
            start_indx(j) = 1;
        end
        Threshold(j) = findoptRth(calcium(start_indx(j):start_indx(j+1),:), Opts);

        %if there is something wrong with finding the threshold
        if isnan(Threshold(j))
            Threshold(j) = mean(Threshold, 'omitnan');
        end

    end

figure
Opts.figs = 0;
Opts.direction_hold = 1;
for j = 1:length(start_indx)-1 %network analysis would be start to start of oscillation
    [N, adj_multi, ~, ~,~,Rij,~] = links(calcium(start_indx(j):start_indx(j+1),:), Threshold(j),Opts,fignum); %%This is where the network is built
    [sorted, cellsor]= sort(N);
    Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);
    degree(:,j) = sum(adj_multi);
    Adj(:,j) = mean(Rij);
    Hubs_multi(:,j) = cellsor(end-top10:end); %Make "hubs" top 10% of highly connected cells
    if i == 1
        out.Hubs_multi(1:length(Hubs),j) = Hubs;
         out.num_hubs(j) = length(Hubs);
         out.degree(:,j) = sum(adj_multi);
         out.top10cells(:,j) = cellsor(end-top10:end);
         out.correlation(j) = mean2(Rij);
    else
         lengthofprevious = size(Hubs_multi,1);
         out.Hubs_multi_pka(1:length(Hubs),j) = Hubs; %so that we can get pka and control
         out.num_hubs(j+lengthofprevious) = length(Hubs);
         out.Adj(:,j+lengthofprevious) = mean(Rij);
         out.degree(:,j+lengthofprevious) = sum(adj_multi);
         out.top10cells_pka(:,j) = cellsor(end-top10:end);
         out.correlation_pka(j) = mean2(Rij);
    end
    %network analysis is performed
  
    %%------plot network------%%
    if 0
    c = 0.7.*ones(size(Locations));
    c(Hubs,1) = ones(length(Hubs), 1).*0;%136/255;
    c(Hubs,2) = ones(length(Hubs), 1).*0;%8./255';
    c(Hubs,3) = ones(length(Hubs), 1).*1;%8./255';
    hold off, plot(graph(sparse(adj_multi)), 'Xdata',Locations(:,1), 'Ydata',Locations(:,2), ...
       'Zdata',Locations(:,3),'NodeColor',c, 'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3, 'EdgeAlpha',0.1)
    %hidden graphs for legend
    hold on,
    ghubs = plot(graph(sparse([0,1;1,0])), 'NodeColor','k',...
        'Xdata',Locations(56:57,1),'Ydata',Locations(56:57,2),'Zdata',Locations(56:57,3),...
        'MarkerSize',15,  'EdgeColor','blue', 'LineWidth',3, 'EdgeAlpha',1);
    hold on,
    gnonhubs = plot(graph(sparse([0,1;1,0])), 'NodeColor',[0.7,0.7,0.7],...
        'Xdata',Locations(56:57,1),'Ydata',Locations(56:57,2),'Zdata',Locations(56:57,3),...
        'MarkerSize',15,  'EdgeColor','k', 'LineWidth',3, 'EdgeAlpha',1);
    legend([ghubs, gnonhubs], {'Hubs','Non-hubs'},'Location', 'northeast')
    set(gca, 'color','none')
    set(gca, 'box','off')
    set(gca,'ytick',[]), set(gca,'xtick',[]), set(gca,'ztick',[])
    %turn off graphs for legend
    set(ghubs,'Visible','off'), set(gnonhubs,'Visible','off')
    if j == 1
        if i == 1
         gif(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'Hubs_WT.gif'])
        else
         gif(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'Hubs_PKA.gif'])
        end
    else
        gif
    end
    end

    close(figure(6))    
    cell_sorted_all(ct,:)= cellsor;
    ct = ct+1;
end
    [loc_mean, loc_spread, Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(length(start_indx)-1, Locations, degree, Adj, Hubs_multi);
    clear Hubs_multi Adj
    
  %% Quantify distance of all cells from center and radially: 
  IsletCenter = mean(Locations);
  dist_from_center = sqrt((Locations(:,1) - IsletCenter(1)).^2 +(Locations(:,2) - IsletCenter(2)).^2 + (Locations(:,3) - IsletCenter(3)).^2);

  % NEED TO QUANTIFY RADIAL: 



    %find top 10% of cells: 
        if i == 1
            out.Threshold_waves = Threshold; %save threshold
            for k = 1:length(start_indx)-1 %calculate the percentage maintained 
                for j = 1:length(start_indx)-1
                    if k < j
                        intravar(k,j) = length(intersect(out.top10cells(:,k),out.top10cells(:,j)))/length(out.top10cells(:,k));
                    else
                        intravar(k,j) = NaN;
                    end
                end
            end
            out.meanLocation = loc_mean;
            out.STLocation = loc_spread;
            out.intravarall = intravar;
            out.intravar = mean2(intravar(~isnan(intravar)));
            out.degree = degree;
        else     
            out.Threshold_waves = Threshold; %save threshold
           for k = 1:length(find(start_indx > cuttime(2)))-1
                for j = 1:length(find(start_indx > cuttime(2)))-1
                    if k < j
                        intravar_pka(k,j) = length(intersect(out.top10cells_pka(:,k), out.top10cells_pka(:,j)))/length(out.top10cells_pka(:,k));
                    else
                        intravar_pka(k,j) = NaN;
                    end
                end
           end
           out.STLocation_pka = loc_spread;
           out.meanLocation_pka = loc_mean;
           out.intravarall_pka = intravar;
           out.intravar_pka = mean2(intravar_pka(~isnan(intravar_pka)));
           out.degree_pka = degree;
            
           for k = 1:length(find(start_indx_hold > cuttime(2)))
                for j = 1:length(find(start_indx_hold < cuttime(2)))
                    intervar(k,j) = length(intersect(out.top10cells(:,j),out.top10cells_pka(:,k)))/length(out.top10cells_pka(:,k));
                end
           end
            out.intervar_all = mean2(intervar);

            out.intervar = mean2(intervar);
        end

        
        out.COG_KL_net(i).trial = COG_KL;
        out.Degree_KL_net(i).trial = Degree_KL;

 %% analyze the trajectory of cells
    [~,cells_sorted] = sort(mean(degree'));
        
   ranking = degree(cells_sorted,:)./(max(degree));
   
   try
       wavenum = length(start_indx)-1;
       numcells = length(adj);
       xvalues = [1:wavenum];
       yvalues = sprintfc('%d',[1:(numcells)]);
       figure, fig1 = heatmap(xvalues, yvalues, ranking)
   catch
       wavenum = length(start_indx);
       numcells = length(adj);
       xvalues = [1:wavenum];
       yvalues = sprintfc('%d',[1:(numcells)]);
       figure, fig1 = heatmap(xvalues, yvalues, ranking)
   end
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   title('Color: Normalized Degree')
   ylabel('Cell Number')
   xlabel('Oscillation Number')




   
   if saveon
       csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'Degree_KL.csv'], Degree_KL)
       csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'COG_KL.csv'], COG_KL)
   end
   clear degree
    if i == 1
        for j = 1:2
         figure(j)
         ax = gca
            text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'Control','FontSize',30)
        end
    else
        for j = 1:2
            figure(j+2)
            ax = gca
                text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'PKa Application','FontSize',30)
        end
    end
    clear degree

end
    calcium_demeaned = (calcium-min(calcium))./(max(calcium)-min(calcium)); 
    time = timestore;
    figure,
     for i = 1:size(cell_sorted_all,1)-1
     nexttile
     plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line2 = plot(time, calcium_demeaned(:,cell_sorted_all(i,1:10)), 'linewidth',1, 'color', 'red')
     hold on, line1 = plot(time, calcium_demeaned(:,cell_sorted_all(i,end-10:end)), 'linewidth',1, 'color', 'blue')

     xline(time(start_indx_hold(i))), xline(time(start_indx_hold(i+1)))
     xlim([time(start_indx_hold(i)), time(start_indx_hold(i+1))])
    title(['Oscillation Number ' num2str(i)])
    axg.data(i) = gca
    legend([line1(1), line2(1)], {'High Degree','Low Degree'})
    set(gca, 'color','none')
    xlabel('Time (s)')
    ylabel('Normalized Ca^{2+} Fluoresence')
    set(gca, 'box','off')

     end
%      saveas(gcf, [savename fileloc(end-5:end-1) filename 'AllOscillations.fig'])
% 

% plot locations
     figure, 
    tt = tiledlayout(2, ceil(size(cell_sorted_all,1)./2))

    tt.TileSpacing = 'compact'
    tt.Padding = 'compact'
        for i =  1:size(cell_sorted_all,1)
            nexttile
            scatter3(Locations(:,1), Locations(:,2), Locations(:,3), 75, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor',[0.7, 0.7, 0.7] , 'MarkerFaceAlpha', 0.5)
            hold on
            scatter3(Locations(cell_sorted_all(i,end-top10:end),1), Locations(cell_sorted_all(i,end-top10:end),2), Locations(cell_sorted_all(i,end-top10:end),3), 100, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue' )

            scatter3(Locations(cell_sorted_all(i,1:top10),1), Locations(cell_sorted_all(i,1:top10),2), Locations(cell_sorted_all(i,1:top10),3), 100, 'MarkerFaceColor', 'red', 'MarkerEdgeColor','red' )
            legend('Normal Cell','High Degree','Low Degree','Location', 'Northeast')
            set(gca, 'color','none')
            title(['Oscillation Number ' num2str(i)])
        
            xlabel('X (\mum)')
            ylabel('Y (\mum)')
            zlabel('Z (\mum)')

        end

        % saveas(gcf, [fileloc 'HighPhaseLocation.png'])
         saveas(gcf, [savename '/' filename 'NetworkLocation.fig'])


try
    saveAllFigsToPPT([savename '/' filename '/Network'])
catch
    mkdir([savename '/' filename])
    saveAllFigsToPPT([savename '/' filename '/Network'])
end
end