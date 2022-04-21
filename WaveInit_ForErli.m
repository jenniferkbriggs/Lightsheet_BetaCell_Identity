%% by Jennifer Briggs 02/28/2022
%This script is modified off of code from Vira and Jenn - calculates phase of cells. 
%Input: calcium wave form
%Output: time of phase difference in ms

%Here we load the file
    calcium = readmatrix('Erli_calcium.csv');%readmatrix('/Volumes/Briggs_2TB/3DIslet/Erli_example.csv'); %change this to be wherever you store your csv
    calcium(1:3,:) = []; %the CSV you have has the first 3 rows as NAN so we remove them
    time = calcium(:,1);   %time is in the first column so pull this out;
    calcium(:,1) = [];     %remove the time so now 'calcium' only has calcium intensity
    figure, plot(time,calcium) %plot calcium

    disp("Resize the figure and then click continue after you are happy with it")
    keyboard
    Locations = readmatrix('Erli_location.csv');
    
%There are a few ways you can look at this wave form - either look at the
%whole thing or look at a single oscillation. Here we define what area we
%are look. Generally, if you are looking for a wave initiator, you'll want
%to look at just the beginning of an oscillation
    title('Select the beginning of 5 oscillations')
    starttime =  ginput(5)  %here you put the time in seconds that you want to start the analysis
    starttime = starttime(:,1);
    hold on, xline(starttime)
    title('Select the peak of those oscillations')
    endtime =    ginput(5)  %put the end time here
    endtime = endtime(:,1);
    xline(starttime, 'label','This is where we start')
    xline(endtime, 'label','This is where we end')


%% Start phase analysis MATLAB works by indexing the datapoints. Therefore, we
    %must find what index the time that you want to look at is.

    for i = 1:5
    start_indx_f = find(abs(time-starttime(i))<0.5); 
    start_indx(i) = start_indx_f(1);
    end_indx_f = find(abs(time-endtime(i))<0.5);
    end_indx(i) = end_indx_f(1);
    end

%% Wave origin analysis
    
    
%in order to make the analysis more accurate, we linearly interpolate between the
%points to artificially increase resolution. We may need to play with this
%as the outputs are not getting very good differentiation between phases.
    calcium_demeaned = (calcium-min(calcium))./(max(calcium)-min(calcium)); 
    %Lets normalize all calcium ranges, assuming the average flouresence value for a cell is a
    %reflection on the staining rather than the actual cell's properties
    cashort = [];
    for i = 1:5
    cashort = calcium_demeaned(start_indx(i):end_indx(i),:);
    step = .0005; %this gives how much to interpolate by
    xq = 1:step:size(cashort,1)+1;
    vq1 = interp1(1:size(cashort,1),cashort,xq);
    timeunits = mean(diff(time))*step*1000; %ms

    
    numcells=size(cashort,2);

    calciumT = (vq1);                           % new calcium time course
    
    % 2. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE
    MeanIslet= (mean(calciumT,2));       % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
    clear vq1 cashort xq camax

    % 3. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET) WITH EACH INDIVIDUAL CELL
    tic
    for j=1:numcells % itterative index for cells
       [c(:,j)]=xcov(calciumT(:,j),MeanIslet,'coeff');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    end
    toc

    [maxCV, maxCL]=max(c);
    clear c

   % 4. PLOTTING SIGNAL, XCOV, AND OUTPUTTING MAX XCOV AND CORRESPONDING TIME LAG
    newmaxCLvec = maxCL-mean(maxCL);
    newmaxCLvec = newmaxCLvec/timeunits; %outputs phase difference in ms


    %this is where you get the final output. phasevecsort gives you the
    %sorted vector of the phase lag compared to the islet mean.
    %cells_sorted is what you are really interested in. This gives the cell
    %index in order from high phase (start earlier) to low
    %phase (start later)
    [phasevecsort(i,:), cells_sorted(i,:)] = sort(newmaxCLvec); 
    
    
    
    end
    
    
     figure, plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
     hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(1:4)), 'linewidth',1, 'color', 'blue')
    hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(end-3:end)), 'linewidth',1, 'color', 'red')
    legend([line1(1), line2(1)], {'High Phase','Low Phase'})

    
    X = Locations(:,1);
    Y = Locations(:,2);
    Z = Locations(:,3);
    figure, scatter3(X,Y,Z,100, newmaxCLvec/max(newmaxCLvec), 'filled')
    colormap hot
    h = colorbar
    set(gca, 'visible', 'off')
    
    %% If you want to watch the calcium oscillations:
    figure, nexttile, 
    
    calcium_av= movmean(calcium_demeaned, 5);
    colormap 
 
    nexttile(2)
    plot(time, calcium_demeaned, 'color',[0.9,0.9,0.9])
    hold on, line1 = plot(time, calcium_demeaned(:,cells_sorted(1:4)), 'linewidth',1, 'color', 'blue')
    hold on, line2 = plot(time, calcium_demeaned(:,cells_sorted(end-3:end)), 'linewidth',1, 'color', 'red')
    xline(1, 'linewidth', 4)

    xwaveinit = X(cells_sorted(1:5));
    ywaveinit = Y(cells_sorted(1:5));
    zwaveinit = Z(cells_sorted(1:5));
    
    my = mean(ywaveinit);
    mx = mean(xwaveinit);
    mz = mean(zwaveinit);

    for i = 1:5:length(calcium);

       nexttile(1)
       scatter3(X,Y,Z,100, calcium_demeaned(i,:), 'filled')
       %annotation('textarrow', mx,my,'String', 'Wave Initiators')
       text(mx+50,my-50,mz, 'Wave Initiators', 'FontSize',20)
       colorbar('south','axislocation', 'in')
       caxis manual
       caxis([0, 1])
       set(gca, 'visible', 'off')

       nexttile(2)
        ax = gca
        axc = ax.Children
        axc(1)
        axc(1).Value = i
        
       drawnow
%        if i == 1
%             keyboard
%         end
 
    end