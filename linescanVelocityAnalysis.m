function linescanVelocityAnalysis(fname, prefs)
%
% line scan analysis - use radon transform on linescan stripes (from RBC 
% shadow) to detect RBCV
% Patrick Drew code from github, adapted by Dori, Kira & Orla, July 2017
% updated May 2018 by Kira 
% NB this code now removes bad frames, so no need to manually delete them
% in ImageJ, if you do manually delete any frames, you must remove them
% from all channels
%
% IMPORTANT: 
% if the line scan was in the wrong direction there is a catch for
% calculating the RBCV, the user must save a .txt file into the folder
% titled 'wrongWay.txt', if this is not there then it will be assumed the
% linescan was the right way. If you do try to extract RBCV from data which
% was scanned the wrongway as you would for rightway data, you may end up
% with no velocity values. 
%
%INPUTS-
%fname = experimental directory(s) containing RBCV.tif file (preprocessing
%line scan images for RBCV, i.e. when line goes through centre of vessel,
%and see RBC shadows as stripes)
%vesselCh = specify which PMT channel the vessel is in, 1 = FITC, 2 = TR,
%the average intensity of the other (calcium) channel will be taken
%NB/ if the user has calcium in their line scan path, they can crop around
%the calcium signal and save this as calcium.tif - if this exists, it will
%be used rather than taking the mean intensity of the non vessel channel,
%and the mean intensity of this cropped section of the image will be used
%instead
%
%OUTPUTS-
%no variables are outputted, the data is extracted in continuous form, and
%saved into the experimental directory, as is a figure
%
%this code requires these functions in your path:
%GetVelocityRadonFig, avgDataOverSlidingWindow, findFolders, pixel4ls, 
%findLocoEvents

%if user doesn't input enough arguments
if nargin<2
    
    %automatically assumes vessel is in channel 2
    prefs.vesselCh = 2; %TR 
    %spike remover:
    %calls a function to clean the trace for RBCV and diam
    prefs.removeSpikes = 1; %1 - will run it through code to remove 
    prefs.plotRemovedSpikes = 1; %1 - will run it through code to remove 
    %ms used to calc the window sz for looking at RBCV
    prefs.windowSz = 40; %ms
    prefs.method = 'mean' %used in RemoveSpikes. 
    %remove bad frames from the calculate velocity trace
    %call on findLocoEvents function to find large jitters in the signal
    %these large noise spikes will be removed 
    %set preferences before call function, as these will be different in
    %this case to when want to detect loco events. 
    prefs.minDist = 100; %min dist between noise spikes
    %whether you want to see the plots for how it has found the noise 
    %will set this to zero, as you will see a plot of the velocity before
    %and after noise removal anyway
    prefs.plotFlag = 1; 
    %for locomotion spikes in the signal that last less than 1 sec are
    %removed, as these are not true loco events, however in this case we
    %want to detect ALL signal spikes to remove them so turn this off 
    prefs.flickerFlag = 0;
    prefs.stdThresh = 1; %for use in clean loco - how many SDs below which is 'rest'
    
end 

%check for existence of the RBCV.tif file
find_tif_file = findFolders(fname, '*RBCV.tif');

%loop through all directories with tif file
for a = 1:size(find_tif_file,2)
    clearvars -except a find_tif_file prefs fname 
    
    %find the local folder containing the tif file
    [expDir,~] = fileparts(find_tif_file{1,a});

     
    %find the other channels
    %check if a calcium.tif file exists 
    %this is the prioritised file type
    find_calc = findFolders(expDir, '*calcium.tif');
    %if there is no calcium.tif file, then take the channel which doesn't
    %have the vessel in it
    if size(find_calc,2) == 0 %if calc doesn't exist
        if prefs.vesselCh == 2
            find_calc=findFolders(expDir, '*ch_1.tif');
        else
            find_calc=findFolders(expDir, '*ch_2.tif');
        end
    else
        %inform user that the calcium.tif file has been detected
        disp('calcium.tif exists...');
    end
    %locomotion is in channel 3
    find_loco = findFolders(expDir, '*ch_3.tif');
    %stimulus is in channel 4
    find_stim = findFolders(expDir, '*ch_4.tif');
    
    %looking for .ini file, which contains parameter info
    find_ini_file = findFolders(expDir, '*.ini');
    x = cell2mat(find_ini_file);
    ini_file = ini2struct(x);
    clear x find_ini_file
        
    %get the number of frames, the width and height of the RBCV ch
    %this has been cropped around the correct part of the scan path, so
    %width and height will be different
    info = imfinfo(find_tif_file{1,a});
    width = info(1).Width;
    height = info(1).Height;
    num_images = numel(info);
    
    %reads the .ini file and extracts the pixel size, fps, and mspline
    %this info is used for creating a time vector
    %call separate function to calculate the pixel size for line scan
    pref.lsNm='*_RBCV.csv';
    [pxsz, linePxs]= pixel4ls(expDir, width, pref); %this is in microns
    fps = str2num(ini_file.x_.frames0x2ep0x2esec); %frames per second (Hz)
    mspline = str2num(ini_file.x_.ms0x2ep0x2eline); %ms per line
    lps = str2num(ini_file.x_.lines0x2eper0x2esecond); %lines per second
    clear inf_file;
    
   
    %take the width and height of the calcium channel
    %may be different as could have been cropped (for calcium.tif)
    info = imfinfo(find_calc{1,1});
    width_calc = info(1).Width;
    height_calc = info(1).Height;
    %take the width and height of loco channel - this can also be used for
    %stim channel, as these have not been cropped at all
    info = imfinfo(find_loco{1,1});
    width_loco = info(1).Width;
    height_loco = info(1).Height;    
    
    %load the tif files for each parameter, across all frames
    %initialise an empty array to store the tiff file in as large dataset
    A = zeros(width,height,num_images);
    calc_ch = zeros(width_calc,height_calc,num_images); 
    loco_ch = zeros(width_loco,height_loco,num_images); 
    stim_ch = zeros(width_loco,height_loco,num_images); 
    
    % Load in all the frames
    % inform user as process can be quite slow
    disp('loading images for all channels...');
    for k = 1:num_images %loop frames
        %disp every 100 frames to show progress
        if (mod(k,100) == 0)
            disp(['Frame ', num2str(k), '/', num2str(num_images)])
        end
        %load the line scan data
        A(:,:,k) = imread(find_tif_file{1,a}, k)';
        %load other channels
        calc_ch(:,:,k) = imread(find_calc{1,1},k)';
        loco_ch(:,:,k) = imread(find_loco{1,1},k)';
        stim_ch(:,:,k) = imread(find_stim{1,1},k)';
    end %end of frames loop 
   
    %reorder the data:
    %put all frames on top of each other, so it becomes a 2D image
    %get the time x space line scan info
    rawLine = reshape(A, [width, height*num_images]);
        
    %also cut the other channels to be same size as linescan RBCV
    %for average traces to keep time/space dims the same
    %take a mean intensity value for each frame
    %stim channel:
    stim_ch_ttt = squeeze(nanmean(nanmean(stim_ch,1),2))';
    %find when stim is on, binarise stim channel:
    stim_ch_ttt = stim_ch_ttt>(mean(stim_ch_ttt)+std(stim_ch_ttt));
    [stim_ch_ttt] = findRogueStims(stim_ch_ttt,fps);
    rawStim = zeros(size(A));
    for b=1:num_images
        %resize the data so the mean intensity for each frame is shaped
        %like the the linescan width and height (per frame)
        rawStim(:,:,b) = repmat(stim_ch_ttt(:,b), [size(A,1), size(A,2)]);
    end
    %then reshape these mean intensity values into the 2D form
    rawStim = reshape(rawStim, [width, height*num_images]);
    clear stim_ch_ttt stim_ch; 
    %repeat reshaping of mean intensity data for calcium channel:
    calc_ch_ttt = squeeze(nanmean(nanmean(calc_ch,1),2))';
    rawCalc = zeros(size(A));
    clear b;
    for b = 1:num_images
        rawCalc(:,:,b) = repmat(calc_ch_ttt(:,b), [size(A,1), size(A,2)]);
    end
    rawCalc = reshape(rawCalc, [width, height*num_images]);
    clear calc_ch_ttt calc_ch;
    %repeat reshaping of mean intensity data for locomotion channel:
    loco_ch_ttt = squeeze(nanmean(nanmean(loco_ch,1),2))';
    rawLoco = zeros(size(A));
    clear b;
    for b = 1:num_images
        rawLoco(:,:,b) = repmat(loco_ch_ttt(:,b), [size(A,1), size(A,2)]);
    end
    rawLoco = reshape(rawLoco, [width, height*num_images]);
    clear loco_ch_ttt loco_ch;

    %explaining the data format:
    
    %width represents the length of the line you drew
    %multiply height x num_images to recreate the line scan image 
        
    %width (i.e. x) is the number of pixels in the line - i.e. the size of 
    %the line you drew
    %height (i.e. y) - each row represents a scan through the line
    %the num_images (i.e. frames) multiplied by length of the y axis is 
    %equivalent to how many times the line was scanned - this is just how 
    %sciscan stores the info, i.e. can put it back into 2 dimensions
        
    %normalise the image so the values are between 0 & 1
    %then every line scan (i.e. across multiple sessions) will be on the 
    %same scale, and comparible 
    rawLine = rawLine-min(min(rawLine));
    rawLine = rawLine/max(max(rawLine));
    
    %save part of the linescan as a figure, so the user can see if this
    %normalising has worked okay
    %take first 100 points as would expect beginning of linescan to be
    %good, as it is how the user decided to attempt a line scan
    figure;
    imagesc(rawLine(:,1:100));
    title('first 100 points of linescan after normalising data'); 
    %save figure as png
    saveas(gcf, fullfile(expDir, 'EG_lsNormalisation.png'));
    close; %close figure
        
    %correction for line scans which were scanned in wrong direction
    %NB may need more complex fix were account for time travelled
    %search for user created text file to indicate linescan is wrong way
    find_ww = findFolders(expDir, 'wrongWay.txt');
    if size(find_ww,2)>=1
        %if the data is the wrong way, then flip it from left to right
        rawLine = fliplr(rawLine);
        %inform the user that their line scan is the wrong way
        disp('line scan wrong way - flipping...');
    end

    % Executes the radon code from
    % https://sites.esm.psu.edu/~pjd17/Drew_Lab/Resources.html 
    % from the paper:
    % Drew PJ, Blinder P, Cauwenberghs G, Shih AY, Kleinfeld D, Rapid 
    % determination of particle velocity from space-time line-scan data 
    % using the Radon transform, Journal of Computational Neuroscience, 
    % 29(1-2):5-11
    
    % NB Drew uses 40ms window but he uses this for arterioles, we needed a
    % larger window for our capillaries to capture more of the streaks
    % this was tested quite extensively, but user can change the window
    % size if it doesn't work for their needs
    % NOTE the window size must be divisible by 4! 
    %NB for HC vs V1 paper - all extracted at 80
    windowsize = round((prefs.windowSz/mspline)/4)*4; %pixels
    
    % This code finds the angle at which the variance is the greatest, 
    % which would be perpendicular to the angle at which the lines actually 
    % are.
    % INPUTS-
    % the line scan data (i.e. 'rawLine'), and the desired window
    % size within which to calculate the RBCV using a radon transform of
    % the stripes within the window 
    % OUTPUTS-
    % thetas - the time varying angle of the space-time image
    % the_t32 - time points of the angle estimates (in lines)
    % spreadmatrix - matix of variances as a function of angles at each 
    % time point
    % transpose the raw line data as function requires time to be in 1st
    % dimension 
    disp('running radon transform function'); 
    [thetasz32,the_tz32,~] = GetVelocityRadon(rawLine',windowsize);
        
    %as the line scan data is averaged over the window size determined and
    %in overlapping steps, the data from the other channels should also be
    %manipulated in the same way
    %i.e. so they are the same size for comparing across channels
    %this function requires time in the 1st dim, so transpose data
    %stim channel:
    [rawStim_tw]=avgDataOverSlidingWindow(rawStim',windowsize);
    rawStim_tw=rawStim_tw';
    %calcium channel:
    [rawCalc_tw]=avgDataOverSlidingWindow(rawCalc',windowsize);
    rawCalc_tw=rawCalc_tw';
    %locomotion channel:
    [rawLoco_tw]=avgDataOverSlidingWindow(rawLoco',windowsize);
    rawLoco_tw=rawLoco_tw';
    
    %calculate seconds per window, useful for timing info
    %used in place of frames per second for line scan data, as scan style
    %is different 
    spw=windowsize*mspline/1000;
    deltat = spw;

    %not totally sure how this bit works but think it takes the 'thetasz32'
    %which is the angle outputs from 'GetVelocityRadon' and once they are
    %'abs' x2 and take away 90, the noise lies at around 25 and less so it
    %identifies where this noise is as a 'thetaEvent'. minDist goes into this 
    %function also, under prefs. This is set to 100 - this is 100 frames
    %not seconds and is normal minDist - i.e. counts events within this
    %time as the same event.
    [thetaEvents] = findLocoEvents(abs(abs(thetasz32)-90)<25,spw,prefs);
    
    %check if any noise was found:
    if ~isempty(thetaEvents)
        %onset of noise spikes
        thetaOnset = thetaEvents(1,:);
        %offset of noise spikes
        thetaOffset = thetaEvents(2,:);
        if ~isempty(find(thetaOnset==0))
            thetaOnset(find(thetaOnset==0))=1;
        end
        %loop the spike events, and replace the noise data with NaNs to remove
        %the bad velocity calculations
        %these are likely in frames were the signal was too poor to see
        %distinct RBC stripes for radon transform
        for ev = 1:size(thetaEvents,2)
            thetasz32(thetaOnset(ev):thetaOffset(ev)+1) = NaN;
        end
    end
    % Calculates the delta t and delta x required to calculate the velocity
    % from the angle.
    % find the length of the line scan in mm:
    deltax = sum(linePxs)/1000;
    % Calculate velocity, this is Eq 4 from the Drew paper.
    vel = (deltax/deltat)*cot(deg2rad(thetasz32));
    %remove values were no info
    vel(vel<0.1)=NaN;
    %filter
    m = nanmean(vel); 
    s = nanstd(vel);
    vel(abs(vel-m) > 3*s) = NaN;
    velocity=vel';
    %check that the velocity is the correct way round, with time in second
    %dimension
    if size(velocity,1)>size(velocity,2)
        velocity=velocity';
    end
    
    %remove spikes from velocity trace if prefs specify
    if prefs.removeSpikes == 1
        disp('Removing spikes from trace');
        velocity = RemoveSpikes(expDir, velocity, prefs);
    end
    
    %flip the time series to way normally prefer, i.e. with time in the
    %last dimension 
    lineangle = thetasz32';
    timepts = the_tz32';
        
    %create time vector for plots
    %in seconds
    time=(the_tz32*mspline)/1000; 
    time=time';
    
    %check that all the data is the same size
    %if not inform the user, and chop the stim, loco and calc so same size
    %as line scan 
    if size(rawStim_tw,2)>size(time,2)
        disp(['data traces are different sizes to line scan channel, ',...
            'cutting to match']); 
        rawStim_tw=rawStim_tw(:,1:size(time,1));
        rawCalc_tw=rawCalc_tw(:,1:size(time,1));
        rawLoco_tw=rawLoco_tw(:,1:size(time,1));
    end
        
    %% name variables same name as in all other codes

    %stim:
    %binarise stim
    stim=rawStim_tw;
    %find when stim is on, binarise stim channel:
    stim = stim > (mean(stim)+std(stim));
    %set first 50 points to zero, in case of initial spike
    stim(1:50) = 0;
    %locomotion:
    %process raw loco data
    movement=rawLoco_tw;
    [locomotion, SDThreshLoco] = cleanLoco(movement, prefs);
    locomotion_raw = abs(diff(movement));
    locomotion_raw = [locomotion_raw,0];
    %calcium:
    calcium=rawCalc_tw;
    
    %make sure correct way round
    if size(the_tz32,1)>size(the_tz32,2)
        the_tz32=the_tz32';
    end
    if size(vel,1)>size(vel,2)
        vel=vel';
    end
    if size(velocity,1)>size(velocity,2)
        velocity=velocity';
    end
    if size(thetasz32,1)>size(thetasz32,2)
        thetasz32=thetasz32';
    end
    if size(lineangle,1)>size(lineangle,2)
        lineangle=lineangle';
    end
    
    %make sure all same size
    if size(the_tz32,2) > size(time,2)
        disp('warning tz32 not same size');
        disp([num2str(size(the_tz32,2)), ' vs ', num2str(size(time,2))]);
        the_tz32=the_tz32(:,1:size(time,2));
    end
    if size(vel,2) > size(time,2)
        disp('warning vel not same size');
        disp([num2str(size(vel,2)), ' vs ', num2str(size(time,2))]);
        vel=vel(:,1:size(time,2));
    end
    if size(velocity,2) > size(time,2)
        disp('warning velocity not same size');
        disp([num2str(size(velocity,2)), ' vs ', num2str(size(time,2))]);
        velocity=velocity(:,1:size(time,2));
    end
    if size(lineangle,2) > size(time,2)
        disp('warning lineangle not same size');
        disp([num2str(size(lineangle,2)), ' vs ', num2str(size(time,2))]);
        lineangle=lineangle(:,1:size(time,2));
    end
    if size(thetasz32,2) > size(time,2)
        disp('warning thetasz32 not same size');
        disp([num2str(size(thetasz32,2)), ' vs ', num2str(size(time,2))]);
        thetasz32=thetasz32(:,1:size(time,2));
    end
    
    disp(['time(s) using lines: ' num2str(time(end))]);
    
    %save variables
    prefs2output = prefs; 
    matfile = fullfile(expDir, 'contData_ls_RBCV');
    save(matfile,'velocity','lineangle','calcium','movement',...
        'locomotion','locomotion_raw','stim','fps','mspline','time',...
        'timepts','rawLine','lps','spw','rawLine', 'pxsz', 'prefs2output', 'SDThreshLoco', '-v7.3');
    
    %% Plot figure
    
    %continuous data traces - linescan image and velocity calculated
    figure;
    screenSz=get(0,'Screensize');
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) screenSz(4)]);
    %line scan image plot
    p1=subplot(5,1,1);
    imagesc(time, [], rawLine); 
    xlim([time(1),time(end)]);
    xlabel('time(s)');
    ylabel('distance along centre of vessel');
    colormap gray
    title('Line Scan Image');
    %stim trials plot
    p2=subplot(5,1,2);
    plot(time, rawStim_tw); 
    xlim([time(1),time(end)]);
    xlabel('time(s)');
    title('Stim Trials');
    %angle calculated from line scan image
    p3=subplot(5,1,3);
    plot(time,thetasz32,'b')% plots the angle at any given time point
    xlim([time(1),time(end)]);
    xlabel('time(s)')
    ylabel('angle, degrees')
    %velocity calculated from angle 
    p4=subplot(5,1,4);
    plot(time,velocity,'r');
    xlabel('time(s)')
    ylabel('velocity, mm/s')
    xlim([time(1),time(end)]);
    p5=subplot(5,1,5);
    plot(time,calcium,'g');
    xlabel('time(s)')
    ylabel('calcium')
    xlim([time(1),time(end)]);
    %link x axis so can zoom into raw image and see velocity etc calculated
    linkaxes([p1,p2,p3,p4,p5],'x');
    %save figure
    saveas(gcf, fullfile(expDir, 'linescan_contTrace_RBCV.png'));
    %close figure
    close;


end %end of looping tif files

end %end of function 
