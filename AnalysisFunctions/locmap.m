function [out] = locmap(pos, phase, per, opts)
%this function is made for 3D analysis of islet data from Erli Jin (Merrins
%Lab) with the goal of understanding the locality of phase initators. 
%Jennifer Briggs 2022

%pos = numcell x 3 array of coordinates
%phase =  ranking of cells by phase--- 1 x numcell array of phase from -1 to 1 for each cell

%per = percentage (0.001 - 0.5) of cells to include in high/low phase
%opts = stucture: figs = 0 for no figures




% Spherical Coordinates: 
[az, el, r] = cart2sph(pos(:,1), pos(:,2), pos(:,3));
%normalize r so the center is 0
%r = r-mean(r)


% high phase: 1) center of gravity location, 2) spread of center of gravity
highphase = find(phase>=1-per);
highphasepos = pos(highphase,:);
highphasecenter = median(highphasepos);


% phase center spread: distance between high phase compared to distance
% between rest of cells
for i = 1:length(pos)
    for j = 1:length(pos)
        distmat(i,j) = sqrt((pos(i,1) - pos(j,1)).^2+(pos(i,2)-pos(j,2)).^2+(pos(i,3)-pos(j,3)).^2);
    end
end

%cell's distance from high phase center
high_dist = sqrt((pos(:,1)-highphasecenter(1)).^2+(pos(:,2)-highphasecenter(2)).^2+(pos(:,3)-highphasecenter(3)).^2);
high_dist = high_dist./max(high_dist);

% Check average radius of high phase cells
high_rad = std(r(highphase))/std(r);
high_az = std(az(highphase))/std(az);
high_el = std(el(highphase))/std(el);

cart_spread = [high_rad, high_az, high_el];
cart_loc = [median(r(highphase)), median(az(highphase)), median(el(highphase))];

out.cart_spread = cart_spread;
out.cart_loc = cart_loc;


%spread of highphasecenter
highphasespread = mean(nonzeros(triu(distmat(highphase, highphase))))/mean(nonzeros(distmat)); %average distance between all cells divided by average distance of all cells from eachother
out.highphasespread = highphasespread;

% low phase: 1) center of gravity location, 2) spread of center of gravity
lowphase = find(phase<=per);
lowphasepos = pos(lowphase,:);
lowphasecenter = median(lowphasepos);

low_dist = sqrt((pos(:,1)-lowphasecenter(1)).^2+(pos(:,2)-lowphasecenter(2)).^2+(pos(:,3)-lowphasecenter(3)).^2);
low_dist = low_dist./max(low_dist);



% distance between high and low phase center of gravity as a function of
% the islet diameter (1 = across whole islet, 0 = right next to each other)

%linear transformation along the slope of the regression between two
%centers of gravity

%center of islet
isletcenter = mean(pos);
[V,D] = eigs([lowphasecenter;isletcenter; highphasecenter]);

pos_new = pos*V(:,1); %gives cell's position along the line of phase difference
pos_new = (pos_new-min(pos_new))/range(pos_new);


%plot
d = jet(length(pos)); %make colormap
[~, ran] = sort(phase); %sort colormap

if opts.figs
    figure, scatter([1:length(pos)], pos_new, 50, d(ran,:), 'filled')
    c = colorbar
    ylabel(c, 'Phase')
    xlabel('Cell Number')
    ylabel('Location along wave')
end

% direction of the wave -> distance from the center of the islet with +
% being begining of wave and - being end of wave. 

out.high_dist = high_dist;
out.low_dist =  low_dist;
out.pos_new = pos_new;
out.highphasecenter = highphasecenter, 
out.lowphasecenter = lowphasecenter
out.V = V
out.D = D

end