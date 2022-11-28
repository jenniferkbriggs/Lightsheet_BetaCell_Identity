function r_mean = oscillationstability(wavenum, Locations, adj_multi, Hubs_multi)
%Assessing the stability of hubs over multiple oscillations

%Questions addressed: 
% 1. What is the spread of hubs in each oscillation?
% 2. How does the 'center of gravity' change?
% 3. How stable is the individual hub?

% ----- Q1 and Q2 ------ %
%Spread of hubs in cartesian   
[az, el, r] = cart2sph(Locations(:,1), Locations(:,2), Locations(:,3));
%normalize radius:
r = (r-min(r))./range(r);
az = (az-min(az))./range(az);
el = (el-min(el))./range(el);
%How much does cartesian locations change with respect to islet
for i = 1:size(Hubs_multi,2)
    az_spread(i) = std(az(~isnan(nonzeros(Hubs_multi(:,i)))));
    el_spread(i) = std(el(~isnan(nonzeros(Hubs_multi(:,i)))));
    r_spread(i) = std(r(~isnan(nonzeros(Hubs_multi(:,i)))));

    az_mean(i) = mean(az(~isnan(nonzeros(Hubs_multi(:,i)))));
    el_mean(i) = mean(el(~isnan(nonzeros(Hubs_multi(:,i)))));
    r_mean(i) = mean(r(~isnan(nonzeros(Hubs_multi(:,i)))));
end

figure,nexttile,
    plot(az_mean,'o-', 'markersize', 10), hold on, plot(el_mean,'o-', 'markersize',10),...
    hold on, plot(r_mean,'o-', 'markersize',10)
legend('\theta','\phi','Radius')
title('Center of Hubs')
ylabel('Normalized coordinate')
nexttile,
    plot(az_spread,'o-', 'markersize', 10), hold on, plot(el_spread,'o-', 'markersize',10),...
    hold on, plot(r_spread,'o-', 'markersize',10)
legend('\theta','\phi','Radius')
title('Spread of Hubs')
xlabel('Oscillation Number')
ylabel('Normalized coordinate')



% ----- Q3 ----- % How often is a hub a hub
hubcell = unique((unique(Hubs_multi)));
hubcell = Hubs_multi(~isnan(Hubs_multi));
for i = 1:length(hubcell)
    perc_time_hub(i) = length(find(Hubs_multi == hubcell(i)))/size(Hubs_multi,2);
end
figure, histogram(perc_time_hub), xlabel('Percent of time a hub cell is a hub cell'), ylabel('Frequency of percent')