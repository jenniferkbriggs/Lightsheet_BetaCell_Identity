%% Plotting Radius
close all
clear all
clc
% Create ROC cuves for different conditions
addpath('~/Documents/GitHub/UniversalCode/')
datadir = '/Volumes/Briggs_10TB/Merrin/Confocal/'
datafiles = dir([datadir '*Retrain_radius.mat']) %location of save from TrainingRadius
datafiles = [datafiles; dir([datadir '*SecondAttempt_radius.mat'])]

datafilespix  = dir([datadir '*RetrainAnalysis.mat'])
datafilespix = [datafilespix; dir([datadir '*SecondAttempt_Analysis.mat'])]

threshold = [.7:.05:.95];

allthresh = .8183; %percent of pixels in the radius that are in correlation treshold before breaking
th = [threshold allthresh]
thpix = 0.8183; %threshold for pixels


%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

allcor = [];
allscore = [];
for kk = 1:length(th)
    figure, gif([datadir 'Threshpics/Radius' 'TH' num2str(th(kk)) '.gif'])
    ct = 1;

    for i = 1:length(datafiles) %loop over islets
         load([datadir datafiles(i).name]) %load in islet
         pixel = load([datadir datafilespix(i).name]);

         isletcor = [];
         isletscore = [];
        for j = 1:length(FinalCordata) %loop over cells
            [TrueCellx TrueCelly] = find(CellMask == j); 
            
            maxedout = 0;
            radius = 1;
            while maxedout == 0
                pixy = circfilled(radius, (NucLoc(j,:)));
                if length(intersect(pixy, [TrueCellx TrueCelly], 'rows')) ~= length(pixy) %then we've hit the maximum radius
                    maxedout = 1;
                else
                    radius = radius +1;
                end
            end
            maxradius = radius;

            [c, indx] = sort(pixel.FinalCordata(j).Correlation);
            scatter(pixel.FinalCordata(j).Pixelsx(indx), pixel.FinalCordata(j).Pixelsy(indx),15, c,'filled'),h = colorbar;
            ylabel(h, 'Correlation with nucleus center')
            gif
           radii = 1;
           score = 1;
           %set threshold
           while score > th(kk)
               %if percent of pixels in radius is greater than threshold, we can increase radius
                [pix_belowth(radii)] = length(find(nonzeros(FinalCordata(j).Correlation(:,radii)) > thpix));
                [tot_pix] = length(nonzeros(FinalCordata(j).Correlation(:,radii)));
                score = pix_belowth(radii)/tot_pix-radii/100;
                radii = radii+1;
           end
            
            est_radius = radii; %only add one because we added one at the 
            %end of the scoring loop (remember we normally add two because 
            %the original analysis that creates the correaltion starts with 
            %index 1 = radius 3

            [estlocations] = circ(est_radius, NucLoc(j,:));
            hold on, scatter(estlocations(:,1), estlocations(:,2), 16, 'k')
            title(['EST: Threshold:' num2str(th(kk)) ' ' datafiles(i).name ': cell ' num2str(j) ])
            
            error(j,i,kk) = maxradius-est_radius;
            trueradius(ct,kk) = maxradius;
            estradius(ct,kk) = est_radius;
           [reallocations] = circ(maxradius, NucLoc(j,:));

            scatter(reallocations(:,1), reallocations(:,2),15, 'r')
            legend('', 'Estimated','Real')
            hold off
            ct = ct+1;
            gif
            %saveas(gcf, [datadir 'Threshpics\' datafiles(i).name 'cell' num2str(j) '.png'])
        end
    end
end

save('FinalRadiusResultsforThr.mat','error','trueradius','estradius')

figure, boxplot((trueradius-estradius),th)
xlabel('Threshold')
ylabel('True Radius - Estimated Radius')
yline(0)
saveas(gcf, [datadir 'FinalErrorResults.fig'])
