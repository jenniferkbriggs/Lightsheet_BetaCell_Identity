%% Roc curves for finding the best correlation treshold
close all
clear all
clc
% Create ROC cuves for different conditions

datadir = '/Volumes/Briggs_10TB/Merrin/Confocal/'
datafiles = dir([datadir '*RetrainAnalysis.mat'])
datafiles = [datafiles; dir([datadir '*SecondAttempt_Analysis.mat'])]

%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

allcor = [];
allscore = [];


for i = 1:length(datafiles) %loop over islets
     load([datadir datafiles(i).name]) %load in islet
     
     isletcor = [];
     isletscore = [];
     
   

    for j = 1:length(FinalCordata) %loop over cells
        [TrueCellx TrueCelly] = find(CellMask == j); 
        
        [c, indx] = sort(FinalCordata(j).Correlation);
%         figure, scatter(FinalCordata(j).Pixelsx(indx), FinalCordata(j).Pixelsy(indx),15, c,'filled'),h = colorbar;
%         ylabel(h, 'Correlation with nucleus center')
%         
%         hold on, scatter(TrueCellx,TrueCelly, 16, 'k')
%         title([datafiles(i).name ': cell ' num2str(j)])
%         saveas(gcf, [datadir datafiles(i).name 'cell' num2str(j) '.fig'])
%         saveas(gcf, [datadir datafiles(i).name 'cell' num2str(j) '.png'])


        %inefficient way to score the matrix
        for p = 1:length(FinalCordata(j).Pixelsx)        
            if isempty(intersect([FinalCordata(j).Pixelsx(p), FinalCordata(j).Pixelsy(p)],[TrueCellx TrueCelly], 'rows'))
                scores(p) = 0;
            else
                scores(p) = 1;
            end
        end
        

       
    [out.islet(i).X(j).data,out.islet(i).Y(j).data,out.islet(i).iT(j).data,out.islet(i).AUC(j).data, out.islet(i).iopth(j).data] = perfcurve(scores,FinalCordata(j).Correlation, 1);
%     figure, plot(out.islet(i).X(j).data, out.islet(i).Y(j).data)
%     title([{[datafiles(i).name ': cell ' num2str(j)]}; {num2str(out.islet(i).AUC(j).data)}])
    %saveas(gcf, [datadir datafiles(i).name ': cell ' num2str(j) 'ROC.png'])
    isletcor = [isletcor; FinalCordata(j).Correlation];
    isletscore = [isletscore; scores'];
    clear scores
    end
    [out.X(i).data,out.Y(i).data,out.T(i).data,out.AUC(i).data, out.opth(i).data] = perfcurve(isletscore, isletcor,1);
    figure, plot(out.X(i).data, out.Y(i).data)
    title([{['Whole Islet:' datafiles(i).name]}; {num2str(out.AUC(i).data)}])
    saveas(gcf, [datadir 'Whole Islet' datafiles(i).name '.png'])

    allcor = [allcor; isletcor];
    allscore = [allscore; isletscore];
end

    [out.all.X,out.all.Y,out.all.T,out.all.AUC, out.all.opth] = perfcurve(allscore, allcor, 1);
    
    figure, plot(out.all.X, out.all.Y)
    title([{'Five Islets'}; {num2str(out.all.AUC)}])
    hold on, 
    
    dist = sqrt((out.all.X).^2+(1-out.all.Y).^2);
    [minny, mini] = min(dist);
    optimalthreshold = out.all.T(mini);%unique(out.all.T(intersect(find(out.all.X == out.all.opth(1)),find(out.all.X == out.all.opth(1)))))
    plot(out.all.X(mini), out.all.Y(mini),'o')
    xlabel('False Positive Rate (pixels counted as part of the cell when they are not)' )
    ylabel('True Positive Rate')
    saveas(gcf, [datadir 'WholeIslet' 'ROC.png'])

    for i = 1:length(datafiles)
        dist = sqrt((out.X(i).data).^2+(1-out.Y(i).data).^2);
        [minny, mini] = min(dist);
        optimalthresholdindislet(i) = out.T(i).data(mini);
        %unique(out.all.T(intersect(find (out.X(i).data == out.opth(i).data(1)),(find(out.Y(i).data...
        %    == out.opth(i).data(2))))))
        
        
        for j = 1:5
            dist = sqrt((out.islet(i).X(j).data).^2+(1-out.islet(i).Y(j).data).^2);
            [minny, mini] = min(dist);
            optimalthresholdindcell(i,j) = out.islet(i).iT(j).data(mini); %unique(out.islet(i).iT(j).data(intersect(find ...
             %   (out.islet(i).X(j).data == out.islet(i).iopth(j).data(1)),(find(out.islet(i).Y(j).data...
             %   == out.islet(i).iopth(j).data(2))))))
        end
    end
    disp(mean([out.all.AUC out.AUC.data]))

    %because I did each islet twice
    incell = [optimalthresholdindcell(1:5,:), optimalthresholdindcell(6:10,:)];
    
    inislet = [optimalthresholdindislet(1:5); optimalthresholdindislet(6:10)];

    