function out = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, calstore, timestore,start_indx, end_indx)
%% Run Network Analysis
% Jennifer Briggs 06/02/2022
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);
if photobleaching
    calcium = detrend_photobleaching(calcium);
end
Threshold(i) = findoptRth(calcium);
end


Threshold = mean(Threshold)
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

if photobleaching
calcium = detrend_photobleaching(calcium);
end
Opts.printasSubplot = 0;
Opts.Subplotnum = 0;
Opts.figs = 1;
Opts.multiseed = 0;
Opts.multiseednum = 1;
fignum = 1

[N, adj, kperc, histArrayPercShort, Rij] = links(calcium,Threshold,Opts,fignum); %%This is where the network is built
[sorted, cellsor]= sort(N);
%network analysis is performed
[L, EGlob, CClosed, ELocClosed, COpen, ELocOpen, nopath]  = graphProperties(adj);
G = graph(adj);
Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);

%%------plot network------%%
    c = 0.7.*ones(size(Locations));
    c(Hubs,1) = ones(length(Hubs), 1).*136/255;
    c(Hubs,2) = ones(length(Hubs), 1).*8./255';
    c(Hubs,3) = ones(length(Hubs), 1).*8./255';
    figure, plot(graph(sparse(adj)), 'Xdata',Locations(:,1), 'Ydata',Locations(:,2), ...
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
       
       
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

if photobleaching
calcium = detrend_photobleaching(calcium);
end
% [~, start_indx] = findpeaks(mean(calcium), 'MinPeakProminence', 0.2)
% xline(time(start_indx))
% start_indx = [1, start_indx]
  if size(calcium,1)<size(calcium,2)
     calcium = calcium';
  end

%%
calcium = (calcium-min(calcium))./range(calcium);
for j = 1:length(start_indx)-1
    if start_indx(j) < 1
        start_indx(j) = 1
    end
    Threshold(j) = findoptRth(calcium(start_indx(j):start_indx(j+1),:), Opts);
end

start_indx(end) = length(calcium);

figure
Opts.figs = 0;
Opts.direction_hold = 1;
for j = 1:length(start_indx)-1 %network analysis would be start to start of oscillation
    [N, adj_multi, ~, ~,~,Rij,~] = links(calcium(start_indx(j):start_indx(j+1),:), quantile(Threshold,0.1),Opts,fignum); %%This is where the network is built
    k(j) = mean2(Rij);
    [sorted, cellsor]= sort(N);
    %network analysis is performed
    Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);
    Hubs_multi(1:length(Hubs),j) = Hubs;
    num_hubs(j) = length(Hubs);
    Adj(:,j) = mean(Rij);
    degree(:,j) = sum(adj_multi);
    %%------plot network------%%
    if 0
    c = 0.7.*ones(size(Locations));
    c(Hubs,1) = ones(length(Hubs), 1).*136/255;
    c(Hubs,2) = ones(length(Hubs), 1).*8./255';
    c(Hubs,3) = ones(length(Hubs), 1).*8./255';
    hold off, plot(graph(sparse(adj_multi)), 'Xdata',Locations(:,1), 'Ydata',Locations(:,2), ...
        'Zdata',Locations(:,3),'NodeColor',c, 'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3, 'EdgeAlpha',0.5)
    %hidden graphs for legend
    hold on,
    ghubs = plot(graph(sparse([0,1;1,0])), 'NodeColor',[136/255,8./255,8./255'],...
        'Xdata',Locations(1:2,1),'Ydata',Locations(1:2,2),'Zdata',Locations(1:2,3),...
        'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3, 'EdgeAlpha',0.5);
    hold on,
    gnonhubs = plot(graph(sparse([0,1;1,0])), 'NodeColor',[0.7,0.7,0.7],...
        'Xdata',Locations(1:2,1),'Ydata',Locations(1:2,2),'Zdata',Locations(1:2,3),...
        'MarkerSize',15, 'EdgeColor','k', 'LineWidth',3, 'EdgeAlpha',0.5);
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
end
    close(figure(6))

%take the minimum number of hubs.  %questionable methodology....
[r_mean,Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(length(start_indx), Locations, degree, Adj, Hubs_multi);

    out.COG_KL_net = COG_KL;
    out.Degree_KL_net = Degree_KL;

 %% analyze the trajectory of cells
    [~,cells_sorted] = sort(mean(degree'));
        
   ranking = degree(cells_sorted,:)./(max(degree));
   
   wavenum = length(start_indx)-1;
   numcells = length(adj);
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);
   figure, fig1 = heatmap(xvalues, yvalues, ranking)
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   title('Color: Normalized Degree')
   ylabel('Cell Number')
   xlabel('Oscillation Number')
   
   csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'Degree_KL.csv'], Degree_KL)
   csvwrite(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'COG_KL.csv'], COG_KL)

   
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
end
saveAllFigsToPPT([savename fileloc(end-5:end-1) filename 'Network'])

end