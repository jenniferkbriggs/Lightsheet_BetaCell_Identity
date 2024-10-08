function [locmean,locspread,Dist_from_cog,COG_KL,Degree_KL] = oscillationstability(wavenum, Locations, degree, Adj, Hubs_multi, phase)
%--Assessing the stability of hubs/highphase over multiple oscillations--%

%Questions addressed: 
% 1. What is the spread of hubs in each oscillation?
% 2. How does the 'center of gravity' change?
% 3. How stable is the individual hub?
% 4. Compare regionality vs. cell
% ----- Q1 and Q2 ------ %
%Spread of hubs in cartesian   
[az, el, r] = cart2sph(Locations(:,1), Locations(:,2), Locations(:,3));
%normalize radius:
r = (r-min(r))./range(r);
az = (az-min(az))./range(az);
el = (el-min(el))./range(el);
%How much does cartesian locations change with respect to islet
wavenum = size(Hubs_multi,2);
Hubs_multi(Hubs_multi == 0) = NaN;
Hubs_nozero = Hubs_multi;
if exist('phase')
for i = 1:size(Hubs_multi,2)
    az_spread(i) = std(az((Hubs_multi(:,i))),'omitnan');
    el_spread(i) = std(el((Hubs_multi(:,i))),'omitnan');
    r_spread(i) = std(r((Hubs_multi(:,i))),'omitnan');

    az_mean(i) = mean(az((Hubs_multi(:,i))),'omitnan');
    el_mean(i) = mean(el((Hubs_multi(:,i))),'omitnan');
    r_mean(i) = mean(r((Hubs_multi(:,i))),'omitnan');
    
    cog(i,:) = median(Locations(Hubs_multi(:,i),:)); %center of gravity - take median in case there are cells scattered a bit. 

 locmean(i,:) = mean(Locations(Hubs_multi(:,i),:));

    cog_std(i,:) = (sum((Locations(Hubs_multi(:,i),:)-cog(i,:)).^2,2))./max((sum((Locations(:,:)-cog(i,:)).^2,2)));
end


%Look at distance traved:
else
    for i = 1:wavenum
    az_spread(i) = std(az(Hubs_multi(:,i)));
    el_spread(i) = std(el(Hubs_multi(:,i)));
    r_spread(i) = std(r(Hubs_multi(:,i)));

    az_mean(i) = mean(az(Hubs_multi(:,i)));
    el_mean(i) = mean(el(Hubs_multi(:,i)));
    r_mean(i) = mean(r(Hubs_multi(:,i)));

    
    cog(i,:) = median(Locations(Hubs_multi(:,i),:)); %center of gravity
    %max distance from cog
 locmean(i,:) = mean(Locations(Hubs_multi(:,i),:));

    cog_std(i,:) = (sum((Locations(Hubs_multi(:,i),:)-cog(i,:)).^2,2))./max((sum((Locations(:,:)-cog(i,:)).^2,2)));
    end
end

locspread = mean(cog_std');


%for plotting ----
    numcells = length(degree);
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
%  ---- ---- ----

figure,nexttile,
    plot(az_mean,'o-', 'markersize', 10), hold on, plot(el_mean,'o-', 'markersize',10),...
    hold on, plot(r_mean,'o-', 'markersize',10)
    legend('\theta','\phi','Radius')
    if exist('phase')
        title('Center Top 5% of Wave Initiators')
    else
        title('Center of Hubs')
    end
    ylabel('Normalized coordinate')
    nexttile,
    plot(az_spread,'o-', 'markersize', 10), hold on, plot(el_spread,'o-', 'markersize',10),...
    hold on, plot(r_spread,'o-', 'markersize',10)
    legend('\theta','\phi','Radius')
    if exist('phase')
        title('Spread Top 5% of Wave Initiators')
    else
    title('Spread of Hubs')
    end
    xlabel('Oscillation Number')
    ylabel('Normalized coordinate')
    
% Regionality - how close is a cell from the center of gravity top X% of hubs?
for i = 1:wavenum
    Dist_from_cog(:,i) = sqrt((Locations(:,1)-cog(i,1)).^2+(Locations(:,2)-cog(i,2)).^2+(Locations(:,3)-cog(i,3)).^2); %cell x wavenum
end
%normalize Dist_from_cog
Dist_from_cog = (Dist_from_cog-(min(Dist_from_cog)))./((max(Dist_from_cog))-(min(Dist_from_cog)));
%Sort distance from center of gravity by average for visualization
[~,Dist_from_cog_idx] = sort(mean(Dist_from_cog,2), 'descend');

figure, 
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);
   heatmap(xvalues, yvalues, Dist_from_cog(Dist_from_cog_idx,:))
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   if exist('phase')
   title('Color: Distance from center of wave intiators')
   else
   title('Color: Distance from center of hubs')
   end
   ylabel('Cell Number')
   xlabel('Oscillation Number')
  

for i = 1:wavenum
    [~,Dist_from_cog_indx(:,i)] = sort((Dist_from_cog(Dist_from_cog_idx,i)),'descend');
end
Dist_from_cog_indx = Dist_from_cog_indx./numcells;

figure, 
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);
   heatmap(xvalues, yvalues, Dist_from_cog_indx)
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   if exist('phase')
      title('Color: Ranked by Distance from center of wave intiators')
   else
        title('Color: Ranked by Distance from center of hubs')
   end
   ylabel('Cell Number')
   xlabel('Oscillation Number')
   
   
% ----- Q3 ----- % How often is a hub a hub
hubcell = unique((unique(Hubs_multi)));
hubcell = Hubs_multi(~isnan(Hubs_multi));
for i = 1:length(hubcell)
    perc_time_hub(i) = length(find(Hubs_multi == hubcell(i)))/size(Hubs_multi,2);
     
end

[~,maximize] = max(size(Adj))
if maximize == 2
    Adj = Adj';
end
norm_degree =  (Adj-(min(Adj)))./((max(Adj))-(min(Adj)));
[~, sort_degree] = sort(mean(Adj,2),'descend')

for i = 1:wavenum
    [~,degree_indx(:,i)] = sort((norm_degree(sort_degree,i)),'descend')
end

degree_indx = degree_indx./numcells;
figure, 
   xvalues = [1:wavenum];
   yvalues = sprintfc('%d',[1:(numcells)]);
   heatmap(xvalues, yvalues, degree_indx)
   ax = gca;
   ax.YDisplayLabels = ydisplayvalues;
   colormap('parula')
   if exist('phase')
       title('Color: Ranked by Cell Phase')
   else
   title('Color: Ranked by Cell Degree')
   end
   ylabel('Cell Number')
   xlabel('Oscillation Number')
   
   
% ----- Q4 ----- %Compare KL divergence
%calculate maximum kl divergence
ct = 1 
while ct < 100
    y = normpdf(degree_indx(:,2));
    x = y(randperm(length(y)));
    average_random(ct) = sum(y.*(log(y)-log(x)));
    ct = ct+1;
end
av_rand = mean(average_random);

    for i = 1:wavenum
        for j = 1:wavenum
            y = normpdf(degree_indx(:,i));
            x = normpdf(degree_indx(:,j));
            Degree_KL(i,j) = sum(y.*(log(y)-log(x)));

            y = normpdf(Dist_from_cog_indx(:,i));
            x = normpdf(Dist_from_cog_indx(:,j));
            COG_KL(i,j) = sum(y.*(log(y)-log(x)));
        end
    end

    Degree_KL = Degree_KL./av_rand; %normalize
    COG_KL = COG_KL./av_rand;
    
    
% ---- Q5 ------ % Does the wave propogate on a specific pole

end


