%% Network Analysis of 3D data
%code for network analysis of 3D data
% Jennifer Briggs 06/02/2022


close all 
clear all
clc
%changing defaults
set(0, 'defaultFigureUnits','normalized', 'defaultFigurePosition', [0.4375 0.1100 0.4675 0.5671]);
set(0,'defaultAxesFontSize',16)




addpath('~/Documents/GitHub/Functional_and_Structural_Networks')

try
addpath('/Users/brigjenn/Documents/GitHub/UniversalCode')
end

fileloc = ('/Volumes/Briggs_10TB/Merrin/Ca_Courses/Singlecelltraces/EJ106/')
filename = 'F08 10G'
fullfile = [fileloc filename]
%Here we load the files
    calcium = readmatrix([fullfile '_Plot.csv']);%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    Locations = readmatrix([fullfile ' pos_Detailed.csv']);
    Locations = Locations(:,1:3);

    figure, plot(mean(calcium'))
    title('Select beginning of second phase, beginning of pka administration, and end of time course')
    [cuttime, ~] = ginput(3)
    close(figure(1))
    calstore = calcium;
    timestore = time;
 fignum = 1
 

for i = 1:2
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

calcium = detrend_photobleaching(calcium);

Threshold(i) = findoptRth(calcium);
end
%% Run Network Analysis
Threshold = mean(Threshold)
for i = 1:2
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

calcium = detrend_photobleaching(calcium);

Opts.printasSubplot = 0;
Opts.Subplotnum = 0;
Opts.figs = 1;
Opts.multiseed = 0;
Opts.multiseednum = 1;

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
hubs_same = length(intersect(Hubs_pka, Hubs_control))/min([length(Hubs_pka), length(Hubs_control)])
hubs_diff = length(setxor(Hubs_pka, Hubs_control))/length(unique([Hubs_pka;Hubs_control]))
c = categorical({'Maintained', 'Different'})
figure, bar(c, [hubs_same, hubs_diff]), ylabel(['Percent of possible hubs maintained or different across conditions'])

c = categorical({'Shortest Path Length', 'Global Efficiency', 'Clustering Coefficient','Local Efficiency','Number of cells with no connections', 'Average Correlation'})
figure, bar(c, [Netinfo_control; Netinfo_pka]), legend('Control','PKa')

saveAllFigsToPPT(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename])


%% analysis for consistency in hubs across oscillations:
for i = 1:2
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);

calcium = detrend_photobleaching(calcium);

figure(7), plot(time,calcium) %plot calcium
disp("Resize the figure and then click continue after you are happy with it")
Opts.direction_hold = 1
Opts.figs = 0;
[start_indx, end_indx] = identify_oscillations(calcium', time, 0)
% [~, start_indx] = findpeaks(mean(calcium), 'MinPeakProminence', 0.2)
% xline(time(start_indx))
% start_indx = [1, start_indx]
calcium = calcium';
close(figure(7))

figure(6),
for j = 1:length(start_indx)-1
    if start_indx(j) < 1
        start_indx(j) = 1
    end
    Threshold(j) = findoptRth(calcium(start_indx(j):start_indx(j+1),:), Opts);
end

figure
for j = 1:length(start_indx)-1 %network analysis would be start to start of oscillation
    [N, adj_multi, ~, ~,~,Rij,~] = links(calcium(start_indx(j):start_indx(j+1),:), median(Threshold),Opts,fignum); %%This is where the network is built
    k(j) = mean2(Rij)
    [sorted, cellsor]= sort(N);
    %network analysis is performed
    Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);
    Hubs_multi(1:length(Hubs),j) = Hubs;
    num_hubs(j) = length(Hubs);
    
    %%------plot network------%%
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
    close(figure(6))

%take the minimum number of hubs.  %questionable methodology....
out = oscillationstability(length(start_indx), Locations, adj_multi, Hubs_multi)
    
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
saveAllFigsToPPT(['/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/3DLightSheet/NetworkAnalysis/' fileloc(end-5:end-1) filename 'OscillationHubs'])

