%% Roc curves
close all
clear all
clc
% Create ROC cuves for different conditions

datadir = '/Volumes/Briggs_10TB/Merrin/Confocal/'
datafiles = dir([datadir '*SecondAttempt_Analysis.mat'])
methodtwo = 0; %finds pixel sensitivity and specificity

%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

thresholds = [0.5:0.05:4];
alldev = [];
allradius = [];
for i = 1:length(datafiles) %loop over islets
     load([datadir datafiles(i).name]) %load in islet
     
     isletdev = [];
     isletradiusscore = [];
    for j = 1:length(FinalCordata) %loop over cells
        [TrueCelly TrueCellx] = find(CellMask == j); %y,x because I flipped the image earlier
        
        [c, indx] = sort(FinalCordata(1).Correlation);
        figure, scatter(FinalCordata(1).Pixelsx(indx), FinalCordata(1).Pixelsy(indx),15, c,'filled'),colorbar;
        hold on, scatter(TrueCellx,TrueCelly, 16, 'k')
        saveas(gcf, [datadir datafiles(i).name 'cell' num2str(j) '.fig'])
        saveas(gcf, [datadir datafiles(i).name 'cell' num2str(j) '.png'])

        %Two ways of assessing accuracy - per pixel or per max radius
        
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
        radiusscore = zeros(size(FinalCordata(j).Correlation,2),1);
        radiusscore(1:maxradius) = 1;
        
        %FinalCordata(j).Correlatin structure is a KxN data file where N is
        %the radius index starting at 3. and K = k are the correlation
        %coefficients between a cell on the circumference and the center
        %filled circle of radius 2
        initialaverage = mean(nonzeros(FinalCordata(j).Correlation(:,1)));
        initialstd = std(nonzeros(FinalCordata(j).Correlation(:,1)));
        for ij = 1:size(FinalCordata(j).Correlation,2)
            deviation(ij) = mean(nonzeros(FinalCordata(j).Correlation(:,ij)))/initialaverage;
            %deviation(ij) = std(nonzeros(FinalCordata(j).Correlation(:,ij)))/initialstd;
        end
        
       if methodtwo  
        allpixels = circfilled(size(FinalCordata(j).Correlation,2)+2, (NucLoc(j,:))) %Because of the image rotation
        figure
        hold on, plot(allpixels(:,1), allpixels(:,2), 'o')
        plot(TrueCellx, TrueCelly, 'o')
        %inefficient way to score the matrix
        for p = 1:length(allpixels)        
            if isempty(intersect(allpixels(p,:),[TrueCellx TrueCelly], 'rows'))
                scores(p) = 0;
            else
                scores(p) = 1;
            end
        end
        
           for th = 1:length(thresholds)
               inradius = find(deviation<thresholds(th));
               if length(inradius) > 1
                    inradius = inradius(1)+2; %true radius is index plus two
               elseif length(inradius) == 0
                   inradius = 0;
               end
               %method 1)
               radius_est(th) = inradius;
               %method 2)
               estpixels = circfilled(inradius, (NucLoc(j,:)));

                   %so ineffiecient
                    for p = 1:length(allpixels)        
                        if isempty(intersect(allpixels(p,:),estpixels, 'rows'))
                            scores_est(p) = 0;
                        else
                            scores_est(p) = 1;
                        end
                    end
                    TN = length(find(scores_est == 0 & scores == 0));
                    FN = length(find(scores_est == 0 & scores == 1));
                    TP = length(find(scores_est == 1 & scores == 1));
                    FP = length(find(scores_est == 1 & scores == 0));

                    if TP + FN ~= length(find(scores==1)) | TN+FP ~= length(find(scores==0)) %sanity check
                        disp('Error in ROC')
                        keyboard
                    end
                    specificity(th) = TN/(FP+TN);%TN/(TN+FN)
                    sensitivity(th) = TP/(TP+FN);
           end
          distfrom1 = sqrt((specificity.^2)+(1-sensitivity).^2)

          optthresh = thresholds(find(distfrom1==min(distfrom1)));
          %sort by x axis
          [k, indx] = unique(specificity)
          for kindx = 1:length(k)
            tempspec = find(specificity == k(kindx));
            [~, maxi] = max(sensitivity(tempspec))
            finalindx(kindx) = tempspec(maxi);
          end


          figure, plot(1-specificity(finalindx), sensitivity(finalindx)), hold on,
          plot(1-specificity, sensitivity, 'o') 

          newX = linspace(0,max(1-specificity(finalindx)),50);
          newY = interp1(1-specificity(finalindx),sensitivity(finalindx),newX);

          AUC = trapz(newX,newY)
          xlabel('1-Specificity')
          ylabel('Sensativity')
          title(num2str(AUC))

       end
    [out.cell(i).X(j).data,out.cell(i).iY(j).data,out.cell(i).iT(j).data,out.cell(i).iAUC(j).data, out.cell(i).iopth(j).data] = perfcurve(radiusscore, deviation,1);
    isletdev = [isletdev; deviation'];
    isletradiusscore = [isletradiusscore; radiusscore];
    end
    [out.X(i).data,out.Y(i).data,out.T(i).data,out.AUC(i).data, out.opth(i).data] = perfcurve(isletradiusscore, isletdev,1);


    alldev = [alldev; isletdev];
    allradius = [allradius; isletradiusscore];
end

    [out.all.X,out.all.Y,out.all.T,out.all.AUC, out.all.opth] = perfcurve(isletradiusscore, isletdev, 1);
