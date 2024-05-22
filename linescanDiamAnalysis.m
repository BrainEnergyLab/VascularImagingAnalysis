function linescanDiamAnalysis(fname,prefs)
%
%function to find the diameter of a vessel from a linescan
%this is less complex then for xyFWHM as the vessel is scanned in a
%straight (vertical) line, and there will be no branches
%will step through in 40ms windows and find mean intensity, then FWHM
%written Dec 2017, by Kira & Orla, updated May 2018 by Kira 
%
%INPUTS-
%fname = directory(s) containing diam.tif - preprocessed vessel in imageJ 
%vesselCh = specify which PMT channel the vessel is in, 1 = FITC, 2 = TR,
%the average intensity of the other (calcium) channel will be taken
%NB/ if the user has calcium in their line scan path, they can crop around
%the calcium signal and save this as calcium.tif - if this exists, it will
%be used rather than taking the mean intensity of the non vessel channel,
%and the mean intensity of this cropped section of the image will be used
%instead
%if no user input, automatically assumes vessel is in channel 2 (TR)
%borderSz = this should be an integer (in pixels) to crop the edges of each
%frame, may be useful if there is noise around the edge (from image
%registration) 
%this is automatically set at 0, unless the user changes this value 
%
%this code requires these functions in your path:
%findFolders, pixel4ls

%if not enough arguments are inputted
if nargin<2
    
    %if need to set pixel size:
    %1 - means you are sending in the pixel size, i.e. you set it on 
    %sciscan and it actually worked
    %0 - means you need to call the magic code to work out the pixel size.
    %as when you set it on sciscan, it dint do anything (i.e. old way)
    prefs.pxFlag = 0; 
    %vessel automatically assumed to be TR
    prefs.vesselCh = 2; %TR 
    %no border for cropping frames
    prefs.borderSz = 0;
    prefs.windowSz = 40; 
    %keep at zero - doesnt work well! 
    prefs.Thresh = 0
    prefs.method = 'mean' %used in RemoveSpikes. 
    
    %calls a function to clean the trace for RBCV and diam1=
    prefs.removeSpikes = 1; %1 - will run it through code to remove 
    prefs.plotRemovedSpikes = 1; %1 - will run it through code to remove  
    prefs.imgThresh = 0.5; %std
    prefs.smoothFactor = 0.05; %how much to smooth by - check this after and alter as needs be
    prefs.plot = 1 %1 to plot the smoothed trace and the original trace.
    prefs.stdThresh = 1; %for use in clean loco - how many SDs below which is 'rest'
end

%check for existence of the diam.tif file
find_tif_file = findFolders(fname, '*diam.tif');

%loop through all directories with tif file
for a=1:size(find_tif_file,2)
   
    clearvars -except a find_tif_file prefs fname 
    
    %find the local folder containing the tif file
    [expDir,~] = fileparts(find_tif_file{1,a});
    
    disp(['Processing: ', num2str(a), '/', num2str(size(find_tif_file,2))]);
    disp(['ExpDir: ', extractAfter(expDir,fname)]);
    
    %find the other channels
    %check if a calcium.tif file exists
    %this is the prioritised file type
    find_calc=findFolders(expDir, '*calcium.tif');
    %if there is no calcium.tif file, then take the channel which doesn't
    %have the vessel in itexpDir
    if size(find_calc,2)==0 %if calc doesn't exist
        if prefs.vesselCh==2
            find_calc=findFolders(expDir, '*ch_1.tif');
        else
            find_calc=findFolders(fname, '*ch_2.tif');
        end
    else
        %inform user that the calcium.tif file has been detected
        disp('calcium.tif exists...');
    end
    %locomotion is in channel 3
    find_loco=findFolders(expDir, '*ch_3.tif');
    %stimulus is in channel 4
    find_stim=findFolders(expDir, '*ch_4.tif');
    
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
    if prefs.pxFlag == 0
        disp('calling magic code to calc px sz, as sciscan dint wrk'); 
        [pxsz, linePxs]= pixel4ls(expDir, width); %this is in microns
    else
        disp('getting px sz from ini file, as sciscan did wrk');
        pxsz = str2num(ini_file.x_.x0x2epixel0x2esz) * 1000000; 
    end
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
    %inform user as process can be quite slow
    disp('loading images for all channels...');
    for k = 1:num_images %loop frames
        %disp every 100 frames to show progress
        if (mod(k,100) == 0)
            disp(['Frame ', num2str(k), '/', num2str(num_images)])
        end
        %load the line scan data
        A(:,:,k) = imread(find_tif_file{1,a}, k)';
        %load other channels
        calc_ch(:,:,k)=imread(find_calc{1,1},k)';
        loco_ch(:,:,k)=imread(find_loco{1,1},k)';
        stim_ch(:,:,k)=imread(find_stim{1,1},k)';
    end %end of frames loop

    %reshape images into continuous time series
    % Put all frames on top of each other, so it becomes a 2D image.
    % get the time x space line scan info
    % add border crop to each frame before reshape 
    if prefs.borderSz>0
        rawIm_ttt = A(prefs.borderSz:width-prefs.borderSz,...
            prefs.borderSz:height-prefs.borderSz,:);
    else
        rawIm_ttt=A;
    end
    rawLine = reshape(rawIm_ttt, [width, height*num_images]);
    % Normalize the image so the values are between 0 & 1
    % Every line scan (i.e. across multiple sessions) will be on the same 
    % scale
    rawLine = rawLine-min(min(rawLine));

    %also cut the other channels to be same size as linescan diam
    %for average traces to keep time/space dims the same
    %take a mean intensity value for each frame
    %stim channel:
    stim_ch_ttt = squeeze(nanmean(nanmean(stim_ch,1),2))';
    %find when stim is on, binarise stim channel:
    stim_ch_ttt = stim_ch_ttt>(mean(stim_ch_ttt)+std(stim_ch_ttt));
    [stim_ch_ttt] = findRogueStims(stim_ch_ttt,fps);
    rawStim = zeros(size(A));
    for b = 1:num_images
        %resize the data so the mean intensity for each frame is shaped
        %like the the linescan width and height (per frame)
        rawStim(:,:,b) = repmat(stim_ch_ttt(:,b), [size(A,1), size(A,2)]);
    end
    %then reshape these mean intensity values into the 2D form
    rawStim = reshape(rawStim, [width, height*num_images]);
    clear stim_ch_ttt stim_ch; 
    %repeat reshaping of mean intensity data for calcium channel:
    calc_ch_ttt=squeeze(nanmean(nanmean(calc_ch,1),2))';
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
        
    %define window size for stepping through linescan
    %match window size to that for velocity
    windowsize = round((prefs.windowSz/mspline)/4)*4; %pixels
    stepsize = round(.25*windowsize); %divide the no of pixels 4
    nlines = size(rawLine,2); %time, i.e. no Lines
    npoints = size(rawLine,1); %space, i.e. Pixels
    nsteps = floor(nlines/stepsize)-3;
        
    %loop through frames/rows of the image to find full width half max
    for k=1:nsteps %step through data 
        
        if (mod(k,1000) == 0)
            disp(['Step ', num2str(k), '/', num2str(nsteps)]);  
        end
        
        %used to create time vector (sampling frequency)
        timepts(k)=1+(k-1)*stepsize+windowsize/2;
        
        % take info from Continuous data
        % loop through data in steps of step size * 4
        % this will overlap by one time box on each loop
        % 0:T/4 time box = stepsize+windowsize/2 - i.e. first box
        % 4 steps later (i.e. 4 time boxes on, 4*stepsize) - 
        % (k-1)*stepsize+windowsize
        clear data_hold_ttt data_hold;
        data_hold_ttt = rawLine(:,1+(k-1)*stepsize:(k-1)*stepsize+windowsize)';

        %call a function to threshold the image for getting diameter
        if prefs.Thresh == 1
            if k == 1
                disp('Thresholding image for FWHM.. slow...');
            end
            [data_hold, ~] = thresholdVesselIm(data_hold_ttt, prefs.imgThresh);
        else
            data_hold = data_hold_ttt;
        end
 
        
        %for the first window - output a figure of the thresholded image so
        %can check it has worked
        if k == 1
            figure;
            subplot(211);
            imagesc(data_hold_ttt);
            title('Raw Data');
            if prefs.Thresh == 1
            subplot(212);
            imagesc(data_hold);
            title('Thresholded Data');
            %save figure
            saveas(gcf, fullfile(expDir, 'linescan_threshImg.png'));
            close;
            end
        end
        
        %extract each row of the linescan, for each frame
        %Smooth the data with the loess method using a span of 50% 
        %of the total number of data points
        clear data_smooth;
%         data_smooth=smooth(nanmean(data_hold,1),0.4,'loess');
        data_smooth = nanmean(data_hold,1)'; 
        %data=smooth(squeeze(rawIm(i,j,:)),0.5,'loess');
        
        % Find the half max value.
        halfMax = (min(data_smooth) + max(data_smooth)) / 2;
        %         %find where the data smooth line (with more points) first crosses the half
        %         %max line and last crosses the line
        %         index1 = find(data_smooth>=halfMax, 1, 'first');
        %         index2 = find(data_smooth>=halfMax, 1, 'last');
        
        %new better way of getting FWHM frames
        %find the nearest pixels to the half way threshold
        [~,I] = sort(abs(halfMax-(data_smooth)),'ascend');
        
        %find which points of FWHM curve to select (i.e. cnt be next to
        %each other) - or wnt get a diam reading 
        for m = 1:size(I,1)
            if abs(I(1) - I(m))>5
                useme = m;
                break
            end
        end
        I = [I(1), I(useme)];
        
        %take the min val as start of curve, and max as end
        index1 = min(I);
        index2 = max(I);
        
        if ~isempty(index1) && ~isempty(index2)
            % FWHM in indexes.
            if prefs.pxFlag == 0
                diameter(k) = sum(linePxs(index1:index2-1));
            else
                diameter(k) = (index2-index1 + 1)*pxsz; %microns
            end
        else
            diameter(k) = NaN;
        end
        
        %suppressing plot outputs
        %if user is debugging and wants to check that it is working, can
        %use the plots below (and turn off the suppression)
        if false
            %output plots every 1000 frames, i.e. to check FWHM is working
            if (mod(k,1000)==0)
                %plot the distribution of the data, 
                %and the indices for fwhm found
                figure;
                subplot(2,1,1);
                plot(data_hold');
                subplot(2,1,2);
                plot(nanmean(data_hold,1),'k');
                hold on;
                plot(data_smooth,'b', 'LineWidth',2);
                plot([index1,index1],[halfMax,halfMax],'go');
                plot([index2,index2],[halfMax,halfMax],'ro');
                title(['window:',num2str(k)]);
            end
        end
        
        %also cut other channels in this window size
        calcium(:,k) = nanmean(nanmean(rawCalc(:,1+(k-1)*stepsize: ...
            (k-1)*stepsize+windowsize),1),2);
        stim(:,k) = nanmean(nanmean(rawStim(:,1+(k-1)*stepsize: ...
            (k-1)*stepsize+windowsize),1),2);
        movement(:,k) = nanmean(nanmean(rawLoco(:,1+(k-1)*stepsize: ...
            (k-1)*stepsize+windowsize),1),2);
        
    end %end of stepping through data


    if prefs.removeSpikes == 1
        disp('Removing spikes from trace');
        diameter = RemoveSpikes(expDir, diameter, prefs);
    end
    
    %call a function to smooth the data (ignoring NaNs)
    [diameter_smooth] = smoothKnit(diameter, prefs); 
    saveas(gcf, fullfile(expDir, 'Smooth.png'));
    close
    

    %create a time vector
    %/1000 to convert from ms to seconds
    time = (timepts*mspline)/1000; 
    %create frame vector     
    frame = [1:size(time,2)];
       
    %inform user how much time is in the final data
    disp(['time(s) using lines: ' num2str(time(end))]);
    
    %process loco and stim
    %stim:
    %binarise stim
    %find when stim is on, binarise stim channel:
    stim = stim > (mean(stim)+std(stim));
    %set first 50 points to zero, in case of initial spike
    stim(1:50) = 0;
    %locomotion:
    [locomotion, SDThreshLoco] = cleanLoco(movement, prefs);
    locomotion_raw = abs(diff(movement));
    locomotion_raw = [locomotion_raw,0];
    
    %save variables
    matfile = fullfile(expDir, 'contData_ls_diam');
    prefs2output = prefs; 
    save(matfile,'time','frame','diameter','diameter_smooth','calcium', ...
        'pxsz','stim','movement','locomotion','locomotion_raw', ...
        'prefs2output','timepts','fps','lps','mspline', 'SDThreshLoco', '-v7.3');
    
    %% plot of continuous data traces
    
    %continuous data traces - linescan image and velocity calculated
    figure;
    screenSz=get(0,'Screensize');
    set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3)/2 ...
        screenSz(4)]);
    %line scan image plot
    p1=subplot(4,1,1);
    plot(time, locomotion_raw); 
    xlim([time(1),time(end)]);
    xlabel('time(s)');
    ylabel('A.U.');
    title('Locomotion'); 
    colormap gray
    title('Line Scan Image');
    %stim trials plot
    p2=subplot(4,1,2);
    plot(time, stim); 
    xlim([time(1),time(end)]);
    xlabel('time(s)');
    title('Stim Trials');
    %velocity calculated from angle 
    p3=subplot(4,1,3);
    plot(time,diameter,'r');
    hold on;
    plot(time,diameter_smooth,'k');
    legend('raw','smooth');
    xlabel('time(s)')
    ylabel('um')
    xlim([time(1),time(end)]);
    p4=subplot(4,1,4);
    plot(time,calcium,'g');
    xlabel('time(s)')
    ylabel('calcium')
    xlim([time(1),time(end)]);
    %link x axis so can zoom into raw image and see velocity etc calculated
    linkaxes([p1,p2,p3,p4],'x');
    %save figure
    saveas(gcf, fullfile(expDir, 'linescan_contTrace_diam.png'));
    %close figure
    close;     
    
end %end of looping through tif files

end %end of function 