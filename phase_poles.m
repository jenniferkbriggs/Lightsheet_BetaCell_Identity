function [rot, rot_pka, rot_intra] = phase_poles(numOscillations, top10, bottom10, top10_pka, bottom10_pka, Locations)

%numOscillations - number of oscillations per trial
%top10           - top x% of cells in first trial
%bottom10           - top x% of cells in first trial
%top10           - top x% of cells in first trial
%top10           - top x% of cells in first trial

%OUTPUTS:
% We output three variables which gives the average normalized distance of change for
% each trial. The maximum distance is 1.
% rot   : - rotation within pre administration
% rot_pka: - rotation within trials pka
% rot_intra - rotation between trials

figure, tt = tiledlayout(2, numOscillations)

for i = 1:(numOscillations)
    %find cells close to the center of gravity
    loctop = Locations(top10(i,:),:);
    locbottom = Locations(bottom10(i,:),:);
    % 
    top_COG = trimmean(loctop, 30); %omit 10% outliers;
    bottom_COG = trimmean(locbottom, 30); 

    %keep cells within some close distance
    keep = find(sqrt((loctop(:,1) - top_COG(:,1)).^2 + (loctop(:,2) - top_COG(:,2)).^2 + (loctop(:,3) - top_COG(:,3)).^2) < 50);
    if length(keep)-2>length(loctop(:,1))
        keyboard
    end
    loctop = loctop(keep,:);

    keep = find(sqrt((locbottom(:,1) - bottom_COG(:,1)).^2 + (locbottom(:,2) - bottom_COG(:,2)).^2 + (locbottom(:,3) - bottom_COG(:,3)).^2)<50);
    locbottom = locbottom(keep,:);


    %find direction of most variance (axis of ca wave)
    [n, p] = pca([locbottom; loctop]); %n is the rotation matrix that rotates the islet so that the axis of the calcium wave is along the 'x' axis. The locations are prenormalized
    if i == 1
        n1 = n;
    end
    toprotated = [loctop*n1];
    bottomrotated = [locbottom*n1];

    Locrotated = Locations*n1;

    
    axis_of_wave(:,i) = n(:,1); 
     ax = gca
    nexttile
    scatter(Locrotated(:,1), Locrotated(:,2))
    hold on
    scatter(toprotated(:,1), toprotated(:,2), 'b', 'filled')
    scatter(bottomrotated(:,1), bottomrotated(:,2), 'r', 'filled')
    if i == 1
    legend('Normal Cell','High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Pre-')
    linkaxes([ax, gca], 'xy')

end

for i = 1:(numOscillations)
    loctop = Locations(top10_pka(i,:),:);
    locbottom = Locations(bottom10_pka(i,:),:);
   % 
   top_COG = trimmean(loctop, 30); %omit 10% outliers;
    bottom_COG = trimmean(locbottom, 30); 

    %keep cells within some close distance
    keep = find(sqrt((loctop(:,1) - top_COG(:,1)).^2 + (loctop(:,2) - top_COG(:,2)).^2 + (loctop(:,3) - top_COG(:,3)).^2) < 50);
    if length(keep)-2 > length(loctop(:,1))
        keyboard
    end
    loctop = loctop(keep,:);
   

    keep = find(sqrt((locbottom(:,1) - bottom_COG(:,1)).^2 + (locbottom(:,2) - bottom_COG(:,2)).^2 + (locbottom(:,3) - bottom_COG(:,3)).^2)<50);
    locbottom = locbottom(keep,:);
    
    %find direction of most variance (axis of ca wave
    [n, p] = pca([loctop; locbottom]);
    rotatedmatrix = [loctop*n; locbottom*n];
    
    axis_of_wave_pka(:,i) = n(:,1);


    % for plotting
        toprotated = [loctop*n1];
    bottomrotated = [locbottom*n1];

    Locrotated = Locations*n1;
    nexttile
     scatter(Locrotated(:,1), Locrotated(:,2))
    hold on
    scatter(toprotated(:,1), toprotated(:,2),'b', 'filled')
    scatter(bottomrotated(:,1), bottomrotated(:,2),'r', 'filled')
    if i == 1
    legend('Normal Cell','High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Post Administration')
        linkaxes([ax, gca], 'xy')

   
end
sgtitle('Rotated according to the first calcium wave axis', 'FontSize',20)


%Since we are in three dimensions, the maximum distance is 3 so we
%normalize:

for i = 1:numOscillations
    for j = 1:numOscillations
    rot(i,j) = sum((axis_of_wave(:,i) - axis_of_wave(:,j)).^2);
    rot_pka(i,j) = sum((axis_of_wave_pka(:,i) - axis_of_wave_pka(:,j)).^2);
    rot_intra(i,j) = sum((axis_of_wave_pka(:,i) - axis_of_wave(:,j)).^2);
    end
end

%find maximum rotation
rot_max = max(max(max([rot, rot_pka, rot_intra])));
for i = 1:50000
topcell = randi(length(Locations),2,1);
bottomcell = randi(length(Locations), 2,1);
[n, p] = pca([Locations(topcell,:); Locations(bottomcell,:)]);

topcell = randi(length(Locations),2,1);
bottomcell = randi(length(Locations), 2,1);
[n2, p] = pca([Locations(topcell,:); Locations(bottomcell,:)]);

rot_max = max(rot_max, sum((n(:,1) - n2(:,1)).^2));
end

if rot_max == max(max(max([rot, rot_pka, rot_intra])))
    disp('Rotation was larger than maximum identified')
end

rot = mean(nonzeros(rot))./rot_max; %zeros are comparing the same oscillation so we remove them
rot_pka = mean(nonzeros(rot_pka))./rot_max; %zeros are comparing the same oscillation so we remove them
rot_intra = mean(mean((rot_intra)))./rot_max; %Don't remove zeros because we don't compare the same oscillations are comparing
end

function plot_rotated2d(numOscillations, top10, bottom10, top10_pka, bottom10_pka, Locations)


figure, tt = tiledlayout(2, numOscillations)
for i = 1:(numOscillations)
    loctop = Locations(top10(i,:),:);
    locbottom = Locations(bottom10(i,:),:);
    

    top_COG = trimmean(loctop, 10); %omit 10% outliers;
    bottom_COG = trimmean(locbottom, 10); 

    %keep cells within some close distance
    keep = find(sqrt((loctop(:,1) - top_COG(:,1)).^2 + (loctop(:,2) - top_COG(:,2)).^2 + (loctop(:,3) - top_COG(:,3)).^2) < 50);
    loctop = loctop(keep,:);

    keep = find(sqrt((locbottom(:,1) - bottom_COG(:,1)).^2 + (locbottom(:,2) - bottom_COG(:,2)).^2 + (locbottom(:,3) - bottom_COG(:,3)).^2)<50);
    locbottom = locbottom(keep,:);

    %find direction of most variance (axis of ca wave)
    
    [n, p] = pca([locbottom; loctop]); %n is the rotation matrix that rotates the islet so that the axis of the calcium wave is along the 'x' axis. The locations are prenormalized
    if i == 1
        n1 = n;
    end
    toprotated = [loctop*n];
    bottomrotated = [locbottom*n];

    Locrotated = Locations*n;

    
    axis_of_wave(:,i) = n(:,1); 
     ax = gca
    nexttile
   %scatter(Locrotated(:,1), Locrotated(:,2))
    hold on
    scatter(toprotated(:,1), toprotated(:,2),'b', 'filled')
    scatter(bottomrotated(:,1), bottomrotated(:,2), 'r', 'filled')
    if i == 1
    legend('High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Pre-')

end

for i = 1:(numOscillations)
    loctop = Locations(top10_pka(i,:),:);
    locbottom = Locations(bottom10_pka(i,:),:);
        
    top_COG = trimmean(loctop, 5); %omit 10% outliers;
    bottom_COG = trimmean(locbottom, 5); 

    %keep cells within some close distance
    keep = find(sqrt((loctop(:,1) - top_COG(:,1)).^2 + (loctop(:,2) - top_COG(:,2)).^2 + (loctop(:,3) - top_COG(:,3)).^2) < 50);
    loctop = loctop(keep,:);

    keep = find(sqrt((locbottom(:,1) - bottom_COG(:,1)).^2 + (locbottom(:,2) - bottom_COG(:,2)).^2 + (locbottom(:,3) - bottom_COG(:,3)).^2)<50);
    locbottom = locbottom(keep,:);

    %find direction of most variance (axis of ca wave
    [n, p] = pca([loctop; locbottom]);
    rotatedmatrix = [loctop*n; locbottom*n];
    
    axis_of_wave_pka(:,i) = n(:,1);


    % for plotting
    toprotated = [loctop*n];
    bottomrotated = [locbottom*n];

    Locrotated = Locations*n;
     nexttile
   % scatter(Locrotated(:,1), Locrotated(:,2))
    hold on
    scatter(toprotated(:,1), toprotated(:,2), 'b', 'filled')
    scatter(bottomrotated(:,1), bottomrotated(:,2),'r', 'filled')
    if i == 1
    legend('High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Post Administration')
end

end
function plot_rotated3d(numOscillations, top10, bottom10, top10_pka, bottom10_pka, Locations)
figure, tt = tiledlayout(2, numOscillations)
for i = 1:(numOscillations)
    loctop = Locations(top10(i,:),:);
    locbottom = Locations(bottom10(i,:),:);
    
    %find direction of most variance (axis of ca wave)
    
    [n, p] = pca([locbottom; loctop]); %n is the rotation matrix that rotates the islet so that the axis of the calcium wave is along the 'x' axis. The locations are prenormalized
    if i == 1
        n1 = n;
    end
    toprotated = [loctop*n];
    bottomrotated = [locbottom*n];

    Locrotated = Locations*n;

    
    axis_of_wave(:,i) = n(:,1); 
     ax = gca
    nexttile
   % scatter3(Locrotated(:,1), Locrotated(:,2), Locrotated(:,3))
    hold on
    scatter3(toprotated(:,1), toprotated(:,2),toprotated(:,3), 'b', 'filled')
    scatter3(bottomrotated(:,1), bottomrotated(:,2),bottomrotated(:,3), 'r', 'filled')
    if i == 1
    legend('Normal Cell', 'High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Pre-')

end

for i = 1:(numOscillations)
    loctop = Locations(top10_pka(i,:),:);
    locbottom = Locations(bottom10_pka(i,:),:);
    
    %find direction of most variance (axis of ca wave
    [n, p] = pca([loctop; locbottom]);
    rotatedmatrix = [loctop*n; locbottom*n];
    
    axis_of_wave_pka(:,i) = n(:,1);


    % for plotting
    toprotated = [loctop*n];
    bottomrotated = [locbottom*n];

    Locrotated = Locations*n;
     nexttile
    %scatter3(Locrotated(:,1), Locrotated(:,2), Locrotated(:,3))
    hold on
    scatter3(toprotated(:,1), toprotated(:,2),toprotated(:,3), 'b', 'filled')
    scatter3(bottomrotated(:,1), bottomrotated(:,2),bottomrotated(:,3), 'r', 'filled')
    if i == 1
    legend('High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Post Administration')
end

end


function plot_normal3D(numOscillations, top10, bottom10, top10_pka, bottom10_pka, Locations)
figure, tt = tiledlayout(2, numOscillations)

for i = 1:(numOscillations)
    loctop = Locations(top10(i,:),:);
    locbottom = Locations(bottom10(i,:),:);
    
    %find direction of most variance (axis of ca wave
    nexttile
    %scatter3(Locations(:,1), Locations(:,2), Locations(:,3))
    hold on
    scatter3(loctop(:,1), loctop(:,2),loctop(:,3), 'b', 'filled')
    scatter3(locbottom(:,1), locbottom(:,2),locbottom(:,3), 'r', 'filled')
    if i == 1
    legend('Normal Cell','High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Pre-')

end

for i = 1:(numOscillations)
    loctop = Locations(top10_pka(i,:),:);
    locbottom = Locations(bottom10_pka(i,:),:);
    
    nexttile
    %scatter3(Locations(:,1), Locations(:,2), Locations(:,3))
    hold on
    scatter3(loctop(:,1), loctop(:,2),loctop(:,3), 'b', 'filled')
    scatter3(locbottom(:,1), locbottom(:,2),locbottom(:,3), 'r', 'filled')
    if i == 1
    legend('Normal Cell','High Phase','Low Phase', 'location','west')
    end
    set(gcf, 'color','white')
    set(gca, 'box','off')
    title('Pre-')
    title('Post Administration')   
end
end


