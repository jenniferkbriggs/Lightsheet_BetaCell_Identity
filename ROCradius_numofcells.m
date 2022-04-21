
close all
clear all
clc

%Find the radial location of training
th = 0.8183; %threshold for pixels

datadir = '/Volumes/Briggs_10TB/Merrin/Confocal/'
datafiles = dir([datadir '*Retrain_radius.mat']) %location of save from TrainingRadius
datafiles = [datafiles; dir([datadir '*SecondAttempt_radius.mat'])]

%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

allcor = [];
allscore = [];

ct = 1;
datafiles = [datafiles; dir([datadir '*SecondAttempt_radius.mat'])]

%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

allcor = [];
allscore = [];

ct = 1;
for i = 1:length(datafiles) %loop over islets
     load([datadir datafiles(i).name]) %load in islet data from TrainingRadius
     
     isletscore = [];
     islettruth = [];
     
    
    for j = 1:length(FinalCordata) %loop over cells
        [TrueCellx TrueCelly] = find(CellMask == j);        
        %find maximum radius:
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
            truth = zeros(1,size(FinalCordata(j).Correlation,2));
            truth(1:maxradius-2) = 1; %set truth radius to 1 for all radii eqal to or under the truth remember radius index starts at 3

       for radii = 1:size(FinalCordata(j).Correlation,2)
            [pix_belowth(radii)] = length(find(nonzeros(FinalCordata(j).Correlation(:,radii)) > th));
            [tot_pix] = length(nonzeros(FinalCordata(j).Correlation(:,radii)));
            radius_score(radii) = (pix_belowth(radii)/tot_pix)-radii/100;

        end
       % allscore = [allscore', radius_score'];

        bestscore(ct) = radius_score(maxradius+2); %radii 1 = radii of 3 so add 2 to get the best score.
        bestpix(ct) = pix_belowth(maxradius+2);
        ct = ct+1;


        

       
    %[out.islet(i).X(j).data,out.islet(i).Y(j).data,out.islet(i).iT(j).data,out.islet(i).AUC(j).data, out.islet(i).iopth(j).data] = perfcurve(truth,allscore, 1);
%    figure, plot(out.islet(i).X(j).data, out.islet(i).Y(j).data)
%    title([{[datafiles(i).name ': cell ' num2str(j)]}; {num2str(out.islet(i).AUC(j).data)}])
    %saveas(gcf, [datadir datafiles(i).name ': cell ' num2str(j) 'ROC.png'])
    isletscore = [isletscore; radius_score'];
    islettruth = [islettruth; truth'];
   % clear allscore
    end
    [out.X(i).data,out.Y(i).data,out.T(i).data,out.AUC(i).data, out.opth(i).data] = perfcurve(islettruth, isletscore,1);
    figure, plot(out.X(i).data, out.Y(i).data)
    title([{['Whole Islet:' datafiles(i).name]}; {num2str(out.AUC(i).data)}])
    saveas(gcf, [datadir 'WholeIslet' datafiles(i).name 'ROC.png'])

    allcor = [allcor; isletscore];
    allscore = [allscore; islettruth];
end

    [out.all.X,out.all.Y,out.all.T,out.all.AUC, out.all.opth] = perfcurve(allscore, allcor, 1);
    
    figure, plot(out.all.X, out.all.Y)
    title([{'Five Islets'}; {num2str(out.all.AUC)}])
    hold on, 
    
    dist = sqrt((out.all.X).^2+(1-out.all.Y).^2);
    [minny, mini] = min(dist);
    optimalthreshold = out.all.T(mini);%unique(out.all.T(intersect(find(out.all.X == out.all.opth(1)),find(out.all.X == out.all.opth(1)))))
    plot(out.all.X(mini), out.all.Y(mini),'o')
    saveas(gcf, [datadir 'AllCells_radiusROC.png'])

    for i = 1:length(datafiles)
        dist = sqrt((out.X(i).data).^2+(1-out.Y(i).data).^2);
        [minny, mini] = min(dist);
        optimalthresholdindislet(i) = out.T(i).data(mini);
        %unique(out.all.T(intersect(find (out.X(i).data == out.opth(i).data(1)),(find(out.Y(i).data...
        %    == out.opth(i).data(2))))))
    end
    
    disp(mean([out.all.AUC out.AUC.data]))

    %because I did each islet twice
    
    inislet = [optimalthresholdindislet(1:5); optimalthresholdindislet(6:10)];

    
    