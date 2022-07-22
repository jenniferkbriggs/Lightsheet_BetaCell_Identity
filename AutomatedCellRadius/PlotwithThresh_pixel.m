%% Roc curves for finding the best correlation treshold
close all
clear all
clc
% Create ROC cuves for different conditions
addpath('~/Documents/GitHub/UniversalCode/')
datadir = '/Volumes/Briggs_10TB/Merrin/Confocal/'
datafiles = dir([datadir '*SecondAttempt_Analysis.mat'])
datafiles = [datafiles; dir([datadir '*RetrainAnalysis.mat'])]
threshold = [ 0.7106    0.9334    0.8844    0.7987    0.7553    0.7082    0.9337    0.9069    0.8064    0.8713];
minthresh = min(threshold);
maxthresh = max(threshold);
allthresh = 0.8183;
th = [minthresh allthresh maxthresh]

allcor = [];
allscore = [];
for kk = 1:3
    figure, gif([datadir 'Threshpics/' 'TH' num2str(th(kk)) '.gif'])

    for i = 1:length(datafiles) %loop over islets
         load([datadir datafiles(i).name]) %load in islet

         isletcor = [];
         isletscore = [];
        for j = 1:length(FinalCordata) %loop over cells
            [TrueCellx TrueCelly] = find(CellMask == j); 

            [c, indx] = sort(FinalCordata(j).Correlation);
            scatter(FinalCordata(j).Pixelsx(indx), FinalCordata(j).Pixelsy(indx),15, c,'filled','linewidth',4),h = colorbar;
            ylabel(h, 'Correlation with nucleus center')
            gif
            est_pix = indx(c>th(kk)); %maybe normalize to correlation through the whole thing

            
            hold on, scatter(FinalCordata(j).Pixelsx(est_pix), FinalCordata(j).Pixelsy(est_pix), 16, 'k')
            title(['EST: Threshold:' num2str(th(kk)) ' ' datafiles(i).name ': cell ' num2str(j) ])
            
%             hold offcc
%             scatter(FinalCordata(j).Pixelsx(indx), FinalCordata(j).Pixelsy(indx),15, c,'filled'),h = colorbar;
%             ylabel(h, 'Correlation with nucleus center')
                
            hold on, scatter(TrueCellx,TrueCelly, 16, 'r', 'filled')
            title(['TRAIN: Threshold:' num2str(th(kk)) ' ' datafiles(i).name ': cell ' num2str(j) ])
            legend('','Estimated','True')
            hold off 
            gif
            %saveas(gcf, [datadir 'Threshpics\' datafiles(i).name 'cell' num2str(j) '.png'])
        end
    end
end
