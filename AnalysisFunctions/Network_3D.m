function out = Network_3D(calcium, cuttime, numtrial, photobleaching, savename, fileloc, filename, Locations, timestore,start_indx, saveon)

%This code calculates all network analysis-related measures: 
% Inputs: Calcium - a txn matrix where t is time and n is cells. 
%         cuttime - indexed times denoting start and stop of each condition (see RunNetworkandWave.m)
%         numtrial - number of conditions - for most of the paper Jin, Briggs et al., this is set to 2 because there was (for example) a control and then a GKa condition
%         photobleaching - binary variable which controls whether to detrend the dataset due to photobleaching (was set to 0 for the manuscript)
%         savename, fileloc, filenmae - paths where the figures are saved
%         Locations - array with x,y,z locations for each cell
%         timestore - just the time array
%         start_indx - indices corresponding to the beginning of each oscillation
%         saveon - whether or not to write csvs (default on)


%% Run Network Analysis for 3D system
% Jennifer Briggs 06/02/2022
ct = 1; %just a counter for oscillations
%options:
if ~exist('saveon') %writing output to csv. If no options then write
    saveon = 1
end

%percent of cells to count as "high degree"
perc=0.1;

%Correct for photobleachine
if photobleaching 
    [~,bo] = detrend_photobleaching(calcium(cuttime(1):cuttime(2),:)); % only remove photobleaching based on control!
    trendline = polyval(bo, [0:length(calcium)-1]');
    calstore = (calcium-trendline)./(trendline-min(min(calcium)));
else
    calstore = calcium;
end

%Identify network threshold
for i = 1:numtrial
calcium = calstore(cuttime(i):cuttime(i+1),:);
time = timestore(cuttime(i):cuttime(i+1),:);
cellnum = min([size(calcium,1), size(calcium,2)]);
Opts.Max = cellnum/8;
Opts.Method = 'Degree'
Opts.avDeg = 7;
Threshold(i) = findoptRth(calcium, Opts);
end

    
%% - Analysis consistency - 
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

%Calculate the network for each oscillation! 
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
    
    %get threshold
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
    %Runs network analysis --- important! - links is found in https://github.com/jenniferkbriggs/Functional_and_Structural_Networks.git

    %Network options
    Opts.direction_hold = 1;
    Opts.printasSubplot = 0;
    Opts.Subplotnum = 0;
    Opts.figs = 1;
    Opts.multiseed = 0;
    Opts.multiseednum = 1;
    fignum = 1;


    [N, adj_multi, ~, ~,~,Rij,~] = links(calcium(start_indx(j):start_indx(j+1),:), Threshold(j),Opts,fignum); %%This is where the network is built


    [sorted, cellsor]= sort(N);
    Hubs = cellsor((sorted - min(sorted))/range(sorted)>.60);
    degree(:,j) = sum(adj_multi);
    Adj(:,j) = mean(Rij);
    Hubs_multi(:,j) = cellsor(end-top10:end); %Make "hubs" top 10% of highly connected cells
    if i == 1
        out.Hubs_multi(1:length(Hubs),j) = Hubs;
         out.num_hubs(j) = length(Hubs);
         deg2 = sum(adj_multi);
         out.degree(:,j) = deg2;
         out.top10cells(:,j) = cellsor(end-top10:end);
         zerodeg = find(degree(:,j) == 0);
         out.bottom10cells(1:length(zerodeg),j) = zerodeg; 
         out.correlation(j) = mean2(Rij);
    else
         lengthofprevious = size(Hubs_multi,1);
         out.Hubs_multi_pka(1:length(Hubs),j) = Hubs; %so that we can get pka and control
         out.num_hubs(j+lengthofprevious) = length(Hubs);
         out.Adj(:,j+lengthofprevious) = mean(Rij);
         deg2 = sum(adj_multi);
         out.degree_pka(:,j) = deg2;
         zerodeg = find(out.degree_pka(:,j) == 0)
         out.top10cells_pka(:,j) = cellsor(end-top10:end);
         out.bottom10cells_pka(1:length(zerodeg),j) = zerodeg;
         out.correlation_pka(j) = mean2(Rij);
    end
  
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
  Radius = max(max(abs([Locations(:,1) - IsletCenter(1), Locations(:,2) - IsletCenter(2), Locations(:,3) - IsletCenter(3)])))
  dist_from_center = sqrt((Locations(:,1) - IsletCenter(1)).^2 +(Locations(:,2) - IsletCenter(2)).^2 + (Locations(:,3) - IsletCenter(3)).^2);
  dist_from_center = dist_from_center./max(dist_from_center); %normalize!
  % NEED TO QUANTIFY RADIAL: 
   
  maxdistance = range(Locations);

    %find top 10% of cells: 
        if i == 1
            out.Threshold_waves = Threshold; %save threshold
            for k = 1:length(start_indx)-1 %calculate the percentage maintained 
                for j = 1:length(start_indx)-1
                    if k < j
                        intravar(k,j) = length(intersect(out.top10cells(:,k),out.top10cells(:,j)))/length(out.top10cells(:,k));
                        out.cog_mov(k,j) = sqrt(sum((loc_mean(k,:) - loc_mean(j,:))./maxdistance).^2);
                        %because of the scale free like nature, there are a
                        %lot of 0 degree cells, so we are looking at all of
                        %the 0 degree cells here! The maximum possible
                        %intersection (=1) is limited by the length of the
                        %shortest non zero array.
                        intravar_bottom(k,j) = length(intersect(out.bottom10cells(:,k),out.bottom10cells(:,j)))/min(length(nonzeros(out.bottom10cells(:,k))), length(nonzeros(out.bottom10cells(:,j))));
                    else
                        intravar(k,j) = NaN;
                        intravar_bottom(k,j) = NaN;
                        out.cog_mov(k,j)  = NaN;
                    end
                end
            end
            %get top10 location: 
            for k = 1:size(out.top10cells,2)
                out.distcenter_top10(:,k) = dist_from_center(out.top10cells(:,k)) 
            end
            out.meanLocation = loc_mean;
            out.STLocation = loc_spread;
            out.intravarall = intravar;
            out.intravar = mean2(intravar(~isnan(intravar)));
            out.intravar_bottom = mean2(intravar_bottom(~isnan(intravar_bottom)));

            out.degree = degree;
            out.dist_from_center = dist_from_center;
        else     
            out.Threshold_waves = Threshold; %save threshold
           for k = 1:length(find(start_indx > cuttime(2)))-1
                for j = 1:length(find(start_indx > cuttime(2)))-1
                    if k < j
                        intravar_pka(k,j) = length(intersect(out.top10cells_pka(:,k), out.top10cells_pka(:,j)))/length(out.top10cells_pka(:,k));
                        intravar_bottom_pka(k,j) = length(intersect(out.bottom10cells_pka(:,k),out.bottom10cells_pka(:,j)))/min(length(nonzeros(out.bottom10cells_pka(:,k))), length(nonzeros(out.bottom10cells_pka(:,j))));
                        out.cog_mov_pka(k,j) = sqrt(sum((loc_mean(k,:) - loc_mean(j,:))./maxdistance).^2);

                    else
                        intravar_pka(k,j) = NaN;
                        intravar_bottom_pka(k,j) = NaN;
                        out.cog_mov_pka(k,j) = NaN;
                    end
                end
           end
           for k = 1:size(out.top10cells_pka,2)
                out.distcenter_top10_pka(:,k) = dist_from_center(out.top10cells_pka(:,k)) 
            end
           out.STLocation_pka = loc_spread;
           out.meanLocation_pka = loc_mean;
           out.intravarall_pka = intravar;
           out.intravar_pka = mean2(intravar_pka(~isnan(intravar_pka)));
           out.intravarbottom_pka = mean2(intravar_bottom_pka(~isnan(intravar_bottom_pka)));
           out.degree_pka = degree;
            out.dist_from_center_pka = dist_from_center;

            
           for k = 1:length(find(start_indx_hold > cuttime(2)))
                for j = 1:length(find(start_indx_hold < cuttime(2)))
                    intervar(k,j) = length(intersect(out.top10cells(:,j),out.top10cells_pka(:,k)))/length(out.top10cells_pka(:,k));
                    intervar_bottom(k,j) = length(intersect(nonzeros(out.bottom10cells(:,j)),nonzeros(out.bottom10cells_pka(:,k))))/min(length(nonzeros(out.bottom10cells_pka(:,k))), length(nonzeros(out.bottom10cells(:,j))));
                end
           end
            out.intervar_all = mean2(intervar);

            out.intervar = mean2(intervar);
            out.intervarbottom = mean2(intervar_bottom);
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
          %  text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'Control','FontSize',30)
        end
    else
        for j = 1:2
            figure(j+2)
            ax = gca
                %text(ax.XLim(end)-range(ax.XLim)/4, ax.YLim(end)-range(ax.YLim)/4,'PKa Application','FontSize',30)
        end
    end

end
    casmooth = smoothdata(calcium); 
    cashort = [];
    calcium_demeaned = (calcium-min(casmooth))./(max(casmooth)-min(casmooth)); 


    time = timestore;
    figure,
     for i = 1:3%size(cell_sorted_all,1)-1
     nexttile
     plot(time, calcium_demeaned, 'color',[0.7,0.7,0.7])
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
     
    saveas(gcf, [savename '/' filename '_NetworkOscillation.png'])
        %Check to make sure number of oscillations is the same
      if find(size(out.degree) ~= size(out.degree_pka))
         error('The number of oscillations in pre and post are not the same')
     end

    %Quantify change in degree between trials
    out.change = mean(abs(diff([mean(out.degree'); mean(out.degree_pka')]./max(max([out.degree, out.degree_pka])))));
 

    %find average top 10% from both trials
    [~,pre_top10] = sort(sum(out.degree'));
    out.top10av = pre_top10(end-top10:end);
    [~,pka_top10] = sort(sum(out.degree_pka'));
    out.top10pka = pka_top10(end-top10:end);
   out.top10percent_avg = length(intersect(out.top10av, out.top10pka))./top10;

% plot locations
     figure, 
    tt = tiledlayout(2, ceil(size(cell_sorted_all,1)./2))

    tt.TileSpacing = 'compact'
    tt.Padding = 'compact'
        for i =  1:size(cell_sorted_all,1)
            nexttile
            scatter3(Locations(:,1), Locations(:,2), Locations(:,3), 75, 'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor',[0.7, 0.7, 0.7] , 'MarkerFaceAlpha', 0.1)
            hold on
            scatter3(Locations(cell_sorted_all(i,end-top10:end),1), Locations(cell_sorted_all(i,end-top10:end),2), Locations(cell_sorted_all(i,end-top10:end),3), 100, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue' )
            scatter3(Locations(cell_sorted_all(i,1:top10),1), Locations(cell_sorted_all(i,1:top10),2), Locations(cell_sorted_all(i,1:top10),3), 100, 'MarkerFaceColor', 'red', 'MarkerEdgeColor','red' )

            legend('Normal Cell','High Degree','Low Degree', 'Location', 'Northeast')


            set(gca, 'color','none')
            title(['Oscillation Number ' num2str(i)])
        
            xlabel('X (\mum)')
            ylabel('Y (\mum)')
            zlabel('Z (\mum)')

        end

         saveas(gcf, [savename '/' filename '_NetworkLocation.fig'])
end