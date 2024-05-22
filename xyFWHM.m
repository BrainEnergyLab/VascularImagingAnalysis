function xyFWHM(fname, vesselCh, prefs, optionsFlag)

%function created by Kira, Dec 2017, updated Dec 2018
%
%exp folder needs to be labelled VR, stim or spont 
%loads in vessel.tif from inputted folder (will error if tif doesn't exist) 
%finds diameter changes using FWHM along vessel (for every frame/branch)
%requires manual user input to draw mask around each of branches 
%(i.e. draw mask around vessel, with little-no background signal within mask
%
%INPUTS-
%fname       : folder name with e.g. vessel.tif file(s) - can be a top dir
%vesselCh    : green (e.g.FITC) = 1, red (e.g.TR) = 2; default 2
%VRflag      : no VR = 0, VR environment on = 1; default 0
%prefs       : structure which contains calc and vessel names
%if user wants to edit, can send in own prefs
%default - prefs.vessel='vessel.tif',prefs.calc='cell_sig.mat'
%optionsFlag : 0 means uses automatic settings, 1 means need to send own
%options in via prefs - NB all options need to be covered, default = 0, see
%further down code for options to use
%
%OUTPUTS- 
%no vars outputted from the function, but it will save some figures and 
%mat files into the experimental directories with tif file in
%will save continuous time series with vessel diameter, loco trace, 
%vis stim, etc., a video to show scan, and a fig to show cont traces

%parameters to specify when call function - vesselCh, VRflag, tifname
%if not specified, default assumptions below
if nargin<2
    vesselCh = 2; %assumes 2 (texas red) by default 
    prefs.vessel = 'vessel.tif'; %default assumes tif file named 'vessel.tif' 
    %default assumes no calc pref, will take mean intensity 
    %NB this can be empty for tif mean intensity, and for separated 3 ch
    %mean intensity
    prefs.calc = []; 
    optionsFlag = 0;
elseif nargin<3
    prefs.vessel = 'vessel.tif'; %default assumes tif file named 'vessel.tif' 
    prefs.calc = []; %default assumes no calc pref, will take mean intensity 
    optionsFlag = 0;
elseif nargin<4
    optionsFlag = 0;
end

if optionsFlag == 0
    %threshold for removing background (std) from raw vessel image (frame1)
    prefs.imgThresh = 0.5; %std
    %set length of skel line to find normal
    %the larger the line, the more accurate the perp line, but
    %the more info lost at either end (limits) of skel
    prefs.skelLineLength = 8; %pixels           
    %sliding time window for averaging FWHM for vessel
    %will use spatial averaging for xy movie
    %define window size for stepping through (in pixels)
    prefs.windowsize=5; %pixels
    %step 1 pixel between time windows (for FWHM averaging)
    prefs.stepsize=2; %pixels
end

%turn off warnings for this function
%or will always warn about polyfit being badly conditioned
warning('off','all')
warning

%call function to find the tif file across the dirs in the fname
find_tif_file = findFolders(fname, prefs.vessel);

%loop through all directories with vessel tif file
for a = 1:size(find_tif_file,2)
    
    %find the tif file from top_dir 
    %save experimental directory as a variable
    [exp_dir,~] = fileparts(find_tif_file{1,a});

    %% load vessel tif images and image info
    
    %looking for .ini file inside experimental directory 
    %.ini file contains parameter info (e.g. size of frame)
    find_ini_file = findFolders(exp_dir, '*.ini');
    x = cell2mat(find_ini_file);
    %convert info into a struct, which contains size info
    ini_file = ini2struct(x);
    %clear old unnecessary variables from workspace
    clear x find_ini_file
    
    %read .ini file and extract useful variables
    pxsz = str2num(ini_file.x_.x0x2epixel0x2esz);  % pixel size (metres)
    pxsz_um = pxsz*1000000;                        % pixel size (um)
    fps = str2num(ini_file.x_.frames0x2ep0x2esec); % frames per second
    %clear old unnecessary variable from workspace
    clear ini_file;
    
    %take tif file name out of a cell and into matrix
    vessel_ch = find_tif_file{1,a};
    

    %load other channel info, calcium (2), locomotion (3), and vis stim (4)
    %load calcium location
    %check prefs to see what form data is stored in
    chFlag = 0;
    if ~isempty(prefs.calc) %load calc info from other source (not tif)
        %load calcium into matrix form
        if strfind(prefs.calc,'csv') %calc data comes from excel
            disp('calcium from excel file');
            calcium = xlsread(fullfile(exp_dir,prefs.calc)); %load excel file
            calcium = calcium(:,2)'; %only take column with calc info in
            %NB/ hard coded as 2, as calc is in 2nd column of excel!!
        elseif strfind(prefs.calc,'mat') 
            disp('calcium from auto ROI detector mat file'); 
            %calc data comes from .mat file
            %depends on what user sends in, i.e. autoROI or cell_sig
            load(fullfile(exp_dir,prefs.calc));
            %if user has sent in pref for cell_sig data, needs to change
            %name to calcium 
            if exist('cell_sig')
                calcium = cell_sig;
            %NB/ this is hard-coded as cell_sort saves calcium as cell_sig!
            end
            clear cell_sig;
        end
    else %if empty, will take mean intensity of tif file from non vessel ch
        if exist(fullfile(exp_dir, 'img_3chs.mat')) < 1 %check for 3 chs
            disp('calcium from mean intensity of entire channel image');
            %if not 3 ch, then load data as norm
            if vesselCh == 1
                calc_ch = cell2mat(findFolders(exp_dir, '*ch_2.tif'));
            elseif vesselCh == 2
                calc_ch = cell2mat(findFolders(exp_dir, '*ch_1.tif'));
            end
        else
            disp(['calcium from mean intensity of entire channel,', ...
            'after separation from vessel']);
            chFlag = 1;
            load(fullfile(exp_dir, 'img_3chs.mat'));
        end
    end
    %locomotion/VR location
    loco_ch = cell2mat(findFolders(exp_dir, '*ch_3.tif'));
    %vis stim location
    stim_ch = cell2mat(findFolders(exp_dir, '*ch_4.tif'));
  
    
    %get image info (from vessel tif file)
    info = imfinfo(vessel_ch);
    width = info(1).Width;    %width of frame (pixels)
    height = info(1).Height;  %height of frame (pixels)
    num_frames = numel(info); %number of frames (frames)
    
    %Initialize an empty array to store the tiff file in
    rawIm = zeros(num_frames,width,height); %vessel channel
    %only need to do for calc if taking raw tif file
    if isempty(prefs.calc)
        calcium = zeros(1,num_frames); %calc channel
    end
    movement = zeros(1,num_frames); %locomotion channel
    stim = zeros(1,num_frames); %vis stim channel
    
    %inform user loading tif files, slow process...
    disp('loading tiff files...');
    if chFlag == 0 %check if has come from 3 channels
        for k = 1:num_frames % loop through frames
            %vessel tif frames
            rawIm(k,:,:) = imread(vessel_ch, k)';
            %load other channels
            %calc (2), loco (3), stim (4)
            if isempty(prefs.calc) %don't load tif file if data in other form
                calcium(:,k) = nanmean(nanmean(imread(calc_ch,k)'));
            end
            
            movement(:,k) = nanmean(nanmean(imread(loco_ch,k)'));
            stim(:,k) = nanmean(nanmean(imread(stim_ch,k)'));
        end
    else
        %hard coded var name as comes from another function
        disp('loading vessel data separated as 3rd channel');
        rawIm = chIm_final{3};
        %hard coded var name as comes from another function
        disp('loading calc data separated from 3 channels');
        calcium_ttt = chIm_final{1};
        for k = 1:num_frames % loop through frames
            movement(:,k) = nanmean(nanmean(imread(loco_ch,k)'));
            stim(:,k) = nanmean(nanmean(imread(stim_ch,k)'));
            calcium(:,k) = nanmean(nanmean(calcium_ttt(k,:,:),2),3); 
        end
    end
    %clear old variables
    clear vessel_ch k calc_ch loco_ch stim_ch calcium_ttt;
    
    %% preprocessing loaded data 
    
    %remove computer calc error from calcium 
    %only do on mean intensity tif or excel (imageJ rois) data
    %already been removed if been through cell sort
    if isempty(prefs.calc) || isempty(strfind(prefs.calc,'csv'))
        calcium=calcium-2^15;
    end
    
    %find when stim is on, binarise stim channel:
    stim=stim>(mean(stim)+std(stim)); 
    %set first 50 points to zero, in case of initial spike
    stim(1:50)=0;
    
    %process loco to extract VR info and raw locomotion
    %auto assume no VR, set VR flag to 0
    VRflag = 0; 
    %if there is VR in the exp dir title, then it will load this and update
    %VR flag to 1
    if regexp(exp_dir, '_VR')
        VRflag = 1; 
        %position shows position within VR
        position = movement;
        %velocity shows speed moved through VR ??
        velocity = diff(movement,1);
        clear movement;
        movement = velocity;
        %filter movement to remove position info
        movement(movement>20) = 0;
        movement(movement<4) = 0;
        %as it has been differentiated, add a 0 to end to make same size as
        %other variables
        movement = [movement, 0];
        velocity = [velocity, 0];

    end
    
    %% check if the vessel is vertical - as it struggles to find the 
    %normal in this case, so will flip 
    
    %check if the vessel is vertical, if so, it will need to be flipped so
    %the code can run
    
    figure;
    title('Check if vessel is vertical...');
    imagesc(squeeze(rawIm(10,:,:)));
    colormap gray;
    
    %user input: check if the vessel is vertical
    button = questdlg('Is the vessel vertical?');
    %yes (1), no (2) or cancel (0) are options - this changes button value
    
    %loop options depending on button user presses (on if happy with skel)
    if strcmp(button,'Cancel')==1 %cancel
        disp('exiting xyFWHM function...'); %notify user exiting function
        return;
    elseif strcmp(button,'Yes')==1 %vessel is vertical - needs flipping
        disp('flipping vessel image...'); 
        rawImttt=zeros(size(rawIm));
        for b = 1:size(rawIm,1)
            rawImttt(b,:,:) = rot90(squeeze(rawIm(b,:,:))); 
        end
        clear rawIm;
        rawIm=rawImttt; clear rawImttt; 
    end
    
    
    %close figure
    close;
    
    %% find skeleton from vessel outline (use 1st vessel tif frame)
    
    %threshold raw image, and remove noise for detecting skeleton
    %just use first image, as all other images are registered to this one
    %in preprocessing
    rawIm_ttt = squeeze(rawIm(10,:,:)); %first frame
    %NB this is hard coded as 0.5, but have found this works well,
    %can remove and send into function via prefs if find need to change
    %input raw vessel image (2D) and preferred threshold (std)
    %will output the new thresholded image, and the automatic skeleton
    [~, skeleton] = thresholdVesselIm(rawIm_ttt, prefs.imgThresh);
    
    %set flag for automatic skeleton at 0
    skelFlag = 0;
    
    %% ask for user input:
    %check if happy with auto skeleton, and how many branches there are...
    
    %fuse the autoskeleton and vessel image (frame 1) together
    skelandvess = imfuse(skeleton,squeeze(rawIm(10,:,:)));
    %display figure with vessel and skeleton 
    %user can check it is acceptable
    figure; title('vessel outline with skeleton');
    imagesc(skelandvess);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %user input: check if happy with skeleton 
    button = questdlg('Happy with skeleton?');
    %yes (1), no (2) or cancel (0) are options - this changes button value
    
    %loop options depending on button user presses (on if happy with skel)
    if strcmp(button,'Cancel')==1 %cancel
        disp('exiting xyFWHM function...'); %notify user exiting function
        return;
    elseif strcmp(button,'No')==1 %not happy with skeleton
        %ask user to draw their own skeleton
        %NB/ it is best to draw this as a closed box which follows 
        %vessel shape
        d = dialog('Position',[300 300 250 150],'Name','Draw Skeleton');
        txt = uicontrol('Parent',d,...
            'Style','text',...
            'Position',[20 80 210 40],...
            'String','Click to draw skeleton');
        btn = uicontrol('Parent',d,...
            'Position',[85 20 70 25],...
            'String','OK',...
            'Callback','delete(gcf)');
        clear d txt btn;
        %clear old (auto) skeleton if user is redoing
        clear skeleton;
        %set new skeleton as the box drawn on the vessel plot by user
        skeleton = roipoly(skelandvess);
        %process this roi to convert to thinned skeleton
        %smooth skel line
        ttt = bwmorph(skeleton,'close',Inf);
        ttt = bwmorph(ttt,'thin',Inf);
        skeleton = bwmorph(ttt, 'bridge');
        %set skelFlag=1 if had to hand draw skeleton
        skelFlag = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %close previous figure opened in this function
    close;
    
    %user can check if new skeleton is acceptable:
    %fuse new skeleton and vessel (frame 1) together
    skelandvess = imfuse(skeleton,squeeze(rawIm(10,:,:)));
    %plot the new skeleton and vessel
    figure; title('vessel outline with final skeleton');
    imagesc(skelandvess);
    %NB/wont have option to exit function now, unless user presses ctrl+c
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % User input: Dialog box to find out how many branches there are 
    % (referring to skeleton to select number of branches)
    
    %strings show text user will see
    prompt = 'Enter number of branches:';
    dlg_title = 'How many branches are there?';
    num_lines = 1;
    %take user input and convert string to number
    answer = inputdlg(prompt,dlg_title,num_lines); 
    num_branches = str2num(answer{1}); %this is num of branches
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % User input: Define where branches are using manual ROI mask selector
    
    %position and title of dialogue box
    d = dialog('Position',[300 300 250 150],'Name','Define vessel regions');
    %position and define inner text
    txt = uicontrol('Parent',d,...
        'Style','text',...
        'Position',[20 80 210 40],...
        'String','Click to draw mask to crop around each branch');
    %position and text of button 
    btn = uicontrol('Parent',d,...
        'Position',[85 20 70 25],...
        'String','OK',...
        'Callback','delete(gcf)');
    clear d txt btn;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear h %clear var before loop (in case already defined)
    %use number of branches specified by user previously 
    %to calculate number of ROI masks needed
    for h = 1:num_branches %loop vessel branches
        
        [Branches(h).masks(:,:),~,~] = roipoly(skelandvess);
        %save the ROI masks created - use later when detecting vessel edges
        
    end
    %close figure
    close;
    
    %extract raw images around each branch (within the roi mask):
    %inform user of loop purpose (as quite slow process)
    disp('Cropping images across all frames...');
    
    clear h i; %clear vars before loop (in case already defined)
    for h = 1:num_branches %loop vessel branches
        
        %disp for each new branch so user knows progress
        disp(['Looping through Branches, Branch ', ...
            num2str(h), '/', num2str(num_branches)]);
        
        %preallocate array with zeros to save memory
        Branches(h).rawIm_mask=zeros(num_frames, size(rawIm,2), size(rawIm,3));
        for i=1:num_frames %loop frames 
            
            %take the rawIm inside the mask for each frame
            clear rawIm_ttt; %clear temp var (in case exists already)
            %create temp var for current frame
            rawIm_ttt=squeeze(rawIm(i,:,:)); 
            %take mean of temp var image - to decide how to threshold
            %background
            rawIm_mean=mean(mean(rawIm_ttt));
            %set areas outside of roi mask, i.e. remove background noise 
            %sometimes works with 0, sometimes with 2^15
            %this means when scan vessel, areas outside of mask wont affect
            %FWHM calculation, as they have been removed 
            if rawIm_mean < 2^15 %if numbers are low, set background to 0
                %take mask for specific branch, and find areas of temp
                %image not within the mask
                rawIm_ttt(~Branches(h).masks) = 0; 
            else
                rawIm_ttt(~Branches(h).masks) = 2^15; %set background 2^15
            end %end of size image mean check
            %save temp image (with background removed) into Branches struct
            Branches(h).rawIm_mask(i,:,:)=rawIm_ttt; 
            
        end %end of frames loop
        
        %also take only the skeleton within the mask, i.e. set rest of
        %frame to zero
        skeleton_ttt=skeleton;
        skeleton_ttt(~Branches(h).masks) = 0; %2^15
        Branches(h).skeleton=skeleton_ttt;
        
        %threshold vessel within mask (for each individual branch)
        %will use this (and individual skel) for rough diameter calculation
        clear rawIm_ttt;
        %create temp variable with only 1st frame 
        %rest of frames are image registered to frame 1
        rawIm_ttt = squeeze(Branches(h).rawIm_mask(1,:,:)); 
        %threshold raw vessel image (frame 1)
        %suppress auto skeleton output, as have already defined skeleton
        [threshIm_ttt, ~] = thresholdVesselIm(rawIm_ttt, prefs.imgThresh);
        
        %get rough diameter of vessel (for each branch):
        %(this will be used for drawing normal line)
        %find max size of rawIm
        maxRawIm = max(size(rawIm_ttt));
        %get a rough area for the vessel
        vesselArea = sum(sum(threshIm_ttt));
        %get length of skeleton
        lengthskel = sum(sum(Branches(h).skeleton));
        %divide area by length to get mean width (pixels)
        roughdiameter = vesselArea/lengthskel;
        %save rough diameter of vessel into structure (will be outputted)
        Branches(h).roughdiam = roughdiameter; %NB this will be overwritten
         clear vesselArea threshIm_ttt;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if lengthskel > 0  %check skeleton exists
            
            clear lengthskel; %no longer need this info
            %find coordinates of skeleton pixels
            %find either end of the skeleton
            clear xend yend; %make sure cleared before allocate values
            [xend, yend] = find(bwmorph(Branches(h).skeleton,'endpoints'));
            %join them together with a contour line
            %use bwtraceboundaries to trace the line from one of endpoints
            %'NW' is the initial search direction (northwest)
            clear contour_skelpts; %make sure cleared before allocate values
            contour_skelpts = bwtraceboundary(Branches(h).skeleton, ...
                [xend(1), yend(1)],'NW');
            
            Branches(h).contour_skelpts = contour_skelpts; 
            
            %num skel pts to loop 
            clear nSkelpts; %make sure cleared before allocate values
            nSkelPts = floor((size(contour_skelpts,1)/2)-1);

            %make sure var name is free before allocate
            clear vid_name writer;
            %specify name of video - individually labelled for branch
            vid_name=['Branch',num2str(h),'_FWHM_vesselScan.avi'];
            %open video writer
            writer = VideoWriter(fullfile(fname,vid_name));
            
            
            %% find branch normal line, and edges of raw image, calc FWHM
            %loop through frames and skel pts,

            %inform user image analysis has started (slow process)
            disp('image analysis for FWHM calc');
            
            for i = 1:num_frames %loop frames
                
                %disp every 100 frames to show progress
                if (mod(i,100) == 0)
                    disp(['Branch ', num2str(h), '/', num2str(num_branches), ...
                        ', Frame ', num2str(i), '/', num2str(num_frames)])
                end
                
                %progress through skeleton pts (minus 1 from start and end)
                %this should be an even number so the pixel of interest is
                %in the middle of the line
                for k = 1:nSkelPts-prefs.skelLineLength %loop skeleton pts
                    %minus # from end as normal line is # in length - so 
                    %would go beyond skeleton if tried to take pts after
                    %this
                    
                    if i == 1 && k == 1 %first iteration of frames and skeleton
                        %only create video for 1st frame (for every branch)
                        %dont need to see every frame scanned
                        open(writer);
                        figure;
                        %proof of concept - user can check skeleton scan
                        %looks okay
                    end
                    
                    %want to take 5 skeleton pts for drawing line through (to 
                    %findnormal) - skeleton pt you want will be the middle pt
                    %skeleton pt you want to look at - this is k (skelpt) +
                    %prefs.skelLineLength/2 (i.e. as is in the middle of the
                    %line)
                    %skelpt=k+(round(prefs.skelLineLength/2));
                    %1st skeleton pt (before pixel of interest) - this is k
                    %to maximise length of line
                    skelb4 = k;
                    %last skel pt (after pixel of interest) - this is k +
                    %prespecified length of skel line to draw norm for
                    skelafter = k+(prefs.skelLineLength-1);
                    
                    %plots to be used if debugging
                    %dont turn these on all the time, will be too many figs
                    if false
                        %plot the part of the skeleton being selected
                        figure; title('selected skeleton pt');
                        imshow(imfuse(skeleton,squeeze(rawIm(i,:,:))));
                        hold on;
                        plot(contour_skelpts(skelb4:skelafter,2),...
                            contour_skelpts(skelb4:skelafter,1),'r');
                    end
                    
                    %to find gradient of line between skelpts, fit line
                    %using 1st order polynomial
                    x = contour_skelpts(skelb4:skelafter,2);
                    y = contour_skelpts(skelb4:skelafter,1);
                    %user can edit poly order here if 1 doesnt work
                    polyOrd = 1;
                    thefit = polyfit(x,y,polyOrd); 
                    
                    %find the min and max values of the line, and use 
                    %polyval to fit the polynomial pts
                    xmin = min(x);
                    xmax = max(x);
                    %create vector of linearly spaced points so can scan 
                    %along vessel length
                    %from min and max of new skeleton (within mask)
                    x1 = linspace(xmin,xmax,10000); 
                    y1 = polyval(thefit,x1);
                    
                    %only plot when debugging, else too many figs!
                    if false
                        
                        %plot original contour line in black, 
                        %with newpolyfit line in red
                        figure; title('orig line, and poly line');
                        plot(x,y,'k');
                        hold on;
                        plot(x1,y1,'r');
                        legend('original','poly'); 
                        
                        %check polyfit line fused over the full skeleton
                        %can see if it is in the right place before find
                        %normal
                        figure; title('skeleton with fitted line');
                        imshow(imfuse(skeleton,squeeze(rawIm(i,:,:))));
                        hold on;
                        plot(x1,y1,'r')
                        
                    end
                    
                    %length of normal line for scanning for FWHM calc
                    %make it twice the size of the rough diam calculated
                    normlength = ceil(roughdiameter)*2; 

                    %catch for if your vessel branch is HORIZONTAL:
                    %and your scan (normal) line should be vertical
                    %as the usual line eqn doesnt work in this case
                    %line eqn loop
                    if range(y1) == 0 %EXACTLY horizontal line
                        
                        %create artifical line coordinates for vertical line
                        %make length of rough diam * 2
                        normx = linspace(x1(round(size(x1,2)/2)), ...
                            x1(round(size(x1,2)/2)),10000); 
                        normy = linspace(y1(round(size(x1,2)/2))- ...
                            normlength,y1(round(size(x1,2)/2))+ ...
                            normlength,10000);
                        
                        %set c (intercept) to zero
                        %check if this value is inf later 
                        c = 0;
                        
                    %catch for if your vessel branch is VERTICAL:
                    %NB as i flip the line it likely won't be vert
                    %and your scan (normal) line should be horizontal
                    %as the usual line eqn doesnt work in this case  
                    elseif range(x1) == 0 %EXACTLY VERTICAL line
                        
                        %take original line, and rotate 90 degrees
                        %and make the norm length specified (2*roughdiam)
                        normy = linspace(y1(2),y1(2),10000);
                        normx = linspace(x1(2)-normlength,x1(2)+normlength,...
                            10000);
                        
                        %set c (intercept) to zero
                        %check if this value is inf later 
                        c=0;
                        
                    else
                        
                        %line equations 
                        %to find perpendicular line to skel polyfit line
                        %calculates the difference between adjacent elements 
                        %of the x and y vectors
                        %so we have a derivative at 2x as many points as we
                        %have skeleton pixels
                        % The normal takes the form y = bx+c
                        %2 is because we want the middle point of the line 
                        %(i.e. skelb4, skelpt, skelafter)
                        
                        %differential to find gradient
                        dy = diff(y1)./diff(x1); 
                        % the slope of the normal 
                        %-1/gradient at that point
                        b = -1/dy(2); 
                        % To find the intercept, c = y-bx.
                        c = y1(round(size(x1,2)/2))-b.*x1(round(size ...
                            (x1,2)/2));
                        
                        %this is the formula for the normal 
                        %use new specified x coordinate to create new y
                        %creates a line that will dissect the skeleton branch
                        %so can measure diameter changes across vessel
                        %by taking the points along the length of the line,
                        %and seeing where they cross the edge of the vessel
                      
                        %make the line the correct size for scanning FWHM
                        normx = linspace(x1(round(size(x1,2)/2))- ...
                            normlength,x1(round(size(x1,2)/2))+ ...
                            normlength,10000);
                        normy = b.*normx + c;
                       
                    end %end of line equation loop 
                    
                    %take the normal line, and plot onto vessel image,
                    %and take corresponding image values
                    %if the intercept is calculated as inf - wont work 
                    %catch to ignore incorrect intercept calculations
                    if ~isinf(abs(c)) 
                        
                        %remove any values outside of limits of image mask 
                        %i.e. particularly any negative values
                        xInd = [find(normx<1),find(normx>maxRawIm)];
                        normx(xInd)=[]; %replace these values with empty
                        normy(xInd)=[];
                        %same for y
                        yInd = [find(normy<1),find(normy>maxRawIm)];
                        normx(yInd)=[]; %replace these values with empty
                        normy(yInd)=[];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %create video to show scan along vessel:
                        %writer was opened previously (in if i==1 && k==1)
                        if i==1
                            %plot of vessel and skeleton, with normal line 
                            %plotted through selected skel pixel pts
                            if k>=1
                                imshow(imfuse(Branches(h).skeleton,...
                                    squeeze(Branches(h).rawIm_mask(i,:,:))));
                                hold on;
                                plot(x1,y1,'r', 'LineWidth',2)
                                plot(normx,normy,'b', 'LineWidth',2);
                                title(['branch:', num2str(h), ', frame:' ...
                                    num2str(i), ', skelpt:' num2str(k)]);
                                frame_ttt = getframe(gcf);
                                frame=frame_ttt.cdata;
                                writeVideo(writer,frame);
                                hold off;
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %coordinates for the normal (intersecting) line
                        x_ttt = round(normx);
                        y_ttt = round(normy);
                        
                        %check the normal line coordinates exist 
                        if ~isnan(nanmean(y_ttt)) 
                            if ~isnan(nanmean(x_ttt))
                                %normal lines are in coordinates, create 
                                %image with zeros where theres no line, 
                                %and ones where there is a line 
                                lineIm_ttt = zeros(size(squeeze(rawIm...
                                    (1,:,:))));
                                
                                for l = 1:size(x_ttt, 2) %loop coords size
                                    %set 1 where normal line lies
                                    lineIm_ttt(y_ttt(l), x_ttt(l)) = 1;
                                end
                                clear l x_ttt y_ttt;
                                
                                %shorten the scan line so it does not
                                %surpass size of mask
                                if size(lineIm_ttt,1)>size(squeeze(rawIm...
                                        (1,:,:)),1)
                                    lineIm_ttt=lineIm_ttt(1:size(squeeze...
                                        (rawIm(1,:,:)),1),:);
                                elseif size(lineIm_ttt,2)>size(squeeze...
                                        (rawIm(1,:,:)),2)
                                    lineIm_ttt=lineIm_ttt(:,1:size(squeeze...
                                        (rawIm(1,:,:)),2));
                                end
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %SWITCH BACK TO USING RAW IMAGE NOW
                                %HAVE THE NORMAL LINE FOR CALCULATING FWHM:
                                %inside raw image (i.e. not thresholded), 
                                %find where normal line == 1
                                %the curve will plot the intensity along 
                                %this line, so you can see the change when
                                %switches from vessel to background
                                %this will be used to find the FWHM
                                Branches(h).data_hold{k,i} = ...
                                    Branches(h).rawIm_mask(i, lineIm_ttt==1);
                                
                            else %if no normal line exists
                                
                                %set to NaN
                                Branches(h).data_hold{k,i} = NaN;
                                
                            end %end of normal line (x coord) NaN check
                        end %end of normal line (y coord) NaN check
                        
                    else %if c (intercept) is NaN - cant find normal line
                        
                        %if normal line is inf, cannot look for edges of vessel
                        %set to NaN
                        Branches(h).data_hold{k,i} = NaN; 
                        
                    end %end of checking if normal line is infinite
                       
                end %end of skel loop
            end %end of frame loop
            
            %close video writer once loop is exited
            close(writer);
            %close video figure
            close;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% spatial averaging, and FWHM calc stepping through windows
            %NB prefs at top used to specify window and step sizes (pixels)
            
            disp('stepping through windows for spatial averaging... slow...'); 

            %calculate number of steps to loop through
            nsteps = floor(nSkelPts/prefs.stepsize)-prefs.windowsize;
            
            %preallocate variables with zeros to save memory & speed up
            vessDiam_ttt = ones(nsteps-prefs.windowsize,num_frames)*NaN;

            %check the raw image (from ROI mask) isn't empty before proceed
            if ~isempty(Branches(h).rawIm_mask)
                
                %inform user its the fwhm calc, as takes a while
                disp(['Branch ', num2str(h)]);
                
                tic; 
                
                %to make sure parfor doesn't error
                delete(gcp('nocreate'));
                cluster = parcluster('local');
                cluster.NumWorkers = 4;
                poolObj = parpool(cluster, cluster.NumWorkers,'IdleTimeout', Inf);
                
                for i = 1:num_frames %loop frames
                    for k = 1:nsteps %loop steps (for smoothing)
                        
                        %disp every 100 frames to show progress
                        if floor(i/100)==(i/100) && k == 1
                            disp(['Branch ', num2str(h), '/', num2str(num_branches),...
                                ', Frame ', num2str(i), '/', num2str(num_frames)])
                        end
                       
                        %ensure that the window fits within the size limits
                        if (k-1)*prefs.stepsize+prefs.windowsize <= size(Branches(h).data_hold,1)
                            
                            %extract the intensity curves for entire window
                            data_ttt = Branches(h).data_hold(1+(k-1)...
                                *prefs.stepsize:(k-1)*prefs.stepsize+...
                                prefs.windowsize,i);
                            
                            %predefine vars to stop error in parfor
                            m = 0; 
                            vessSz_um = zeros(1,size(data_ttt,1)); 
                            parfor m = 1:size(data_ttt,1)
                                %refer to the variable in it's unsliced
                                %form to force it to be "broadcast"
                                % https://uk.mathworks.com/matlabcentral/answers/224150-index-exceeds-matrix-dimensions-error-when-using-parfor
                                if m <= size(data_ttt,1) %refer to var not sliced
                                    %smooth data
                                    data_smooth = smooth(data_ttt{m,:},...
                                        0.13,'loess')';
                                    % Find the half max value.
                                    halfMax = (min(data_smooth) + ...
                                        max(data_smooth)) / 2;
                                    
                                    %find where the data smooth line (with more points) first crosses the half
                                    %max line and last crosses the line
                                    index1 = find(data_smooth>=halfMax, 1, 'first');
                                    index2 = find(data_smooth>=halfMax, 1, 'last');
                                    
                                    %check the curve existed so can find dist
                                    if ~isempty(index1) && ~isempty(index2)
                                        %multiply by pixel size to get in um
                                        vessSz_um(m) = (index2 - index1)* pxsz_um;
                                    else
                                        vessSz_um(m) = NaN;
                                    end
                                end %end of forcing code to not slice var
                            end %end of parallel for loop
                            
                            %save into fwhm pixel variable - for every frame, 
                            %and every skel pt
                            vessDiam_ttt(k,i) = nanmean(vessSz_um);
                            clear vessSz_um; 

                        end %end of check if enough space for skel step
                    end %end of looping through steps 
                end %end of looping through frames
                
                toc; 
                delete(poolObj); 
                
            else
                
                %notify user that the vessel image is empty - bad signal
                disp(['Branch ', num2str(h), ': bad signal, not processed']);
                
            end %end of checking if vessel image (from mask) is empty
            
            Branches(h).fwhm_um = vessDiam_ttt;
            clear vessDiam_ttt; 

            %could check for outlier diameters and replace with NaNs? 
            %do for each individual branch (as branches could be diff szs),
            %also do within a frame? as vess dilating may come out as an
            %outlier
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% cont plot - will save into exp folder so user can eyeball data
            
             
            %create frame vector and time vector (for plotting)
            %will be outputted when save
            frames = 1:num_frames;
            time = frames/fps;
            
            %find size of screen (for figure sizing)
            screenSz=get(0,'Screensize');
            
            %plot continuous variables
            figure;
            %make fig size of screen
            set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
                screenSz(4)]);
            %plot vessel diam changes (pixels)
            subplot(511);
            imagesc(frames, [], Branches(h).fwhm_um);
            xlim([frames(1) frames(end)]);
            title(['Continuous Diameter Changes, Branch#',num2str(h)]);
            ylabel('Skel#'); xlabel('Frames');
            %plot average diam change (pixels)
            subplot(512);
            plot(frames, nanmean(Branches(h).fwhm_um,1),'r');
            xlim([frames(1) frames(end)]);
            title('Continuous AVERAGE Diameter Changes');
            xlabel('Frames'); ylabel('um');
            %vis stim 
            subplot(513);
            plot(frames, stim, 'k');
            xlim([frames(1) frames(end)]);
            title('Stim Trials');
            xlabel('Frames');
            %locomotion
            subplot(514);
            plot(frames, movement, 'b');
            xlim([frames(1) frames(end)]);
            title('Continuous Locomotion');
            xlabel('Frames');
            %calcium 
            subplot(515);
            plot(frames, nanmean(calcium,1),'g');
            xlim([frames(1) frames(end)]);
            title(['Continuous Average Calcium, nROIs=', ...
                num2str(size(calcium,1))]);
            xlabel('Frames');
            %save figure as png
            saveas(gcf, fullfile(exp_dir, ['Branch',num2str(h),...
                '_contPlots_fwhm']), 'png');
            close; %close figure
            
            %plot branch onto raw image
            img_ttt=imfuse(Branches(h).skeleton,...
                squeeze(Branches(h).rawIm_mask(1,:,:)));
            figure;
            imagesc(imfuse(img_ttt, squeeze(rawIm(1,:,:))));
            title(['Branch',num2str(h),' overlaid onto raw img']);
            saveas(gcf, fullfile(exp_dir, ...
                ['Branch',num2str(h),'_imagePlot']), 'png');
            close; %close figure
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else
            
            %will skip a branch if no skel created
            disp(['Branch ', num2str(h), ' has no skeleton...']);
               
        end %end of check for existence of branch skeleton
        
    end %end of branch loop
    
    %display avg vess size
    for z = 1:num_branches %loop vessel branches
        disp(['Branch ', num2str(z), ' = ', num2str(nanmean(Branches(z).fwhm_um(:)))]); 
        Branches(z).roughdiam = nanmean(Branches(z).fwhm_um(:)); 
    end
    
    %% make cont trace across branches
    
    if num_branches > 1 %check if more than one branch
        
        %make the xy diam data into cont trace by calling function
        disp('making diam data cont trace...');
        %     [cont_diam] = xyDiamStruct2contTrace(exp_dir);
        
        for b = 1:size(Branches,2) %loop branches
            %find the index of when the entire branches swap (not
            %necessarily most responsive pt of branch)
            if b == 1
                %index of when the data switches from one branch to next
                %in event data is stored continuously (not separated in
                %struct) in case of first branch, just need size of branch1
                branchswapInd(b,:) = [1, size(Branches(b).fwhm_um,1)];
            else
                %in case of branches beyond 1st, find start position by
                %size of previous branch, and end position by adding this
                %to size of current branch
                branchswapInd(b,:) = [branchswapInd(b-1,2)+1, ...
                    (branchswapInd(b-1,2))+(size(Branches(b).fwhm_um,1))];
            end
        end
        
        start = 2; %start of vessel skeleton
        stop  = branchswapInd(end)-1; %end of vessel skeleton
        
        for b = 1:size(Branches,2) %loop branches
            if b == 1 %check if you are on the first branch
                %check if the end of the index is greater than the start pt
                %i.e. is this branch included as responsive
                if start <= branchswapInd(b,2) && ...
                        stop >= branchswapInd(b,1)
                    %check if the index is less than the stop point, i.e.
                    %this whole branch is included
                    %the index is less than the stop pt, so it may continue
                    %to next branch - just take until end of switch pt
                    if branchswapInd(b,2) <= stop
                        cont_diam = Branches(b).fwhm_um ...
                            (start:branchswapInd(b,2),:);
                        %the index is larger than the stop pt, so the stop
                        %pt is within this branch - just use stop pt
                    else
                        cont_diam = Branches(b).fwhm_um(start:stop,:);
                    end
                else %first branch isn't included
                    cont_diam = [];
                end
            else %subsequent branches
                if b >= 2 %branch 2
                    %start is within index, but stop pt goes beyond, so
                    %just take end of branch
                    if start >= branchswapInd(b,1) && ...
                            stop >= branchswapInd(b,2)
                        cont_diam = [cont_diam; Branches(b). ...
                            fwhm_um(start-branchswapInd(b,1):end,:)];
                        %start from beginning of branchswapInd (as start pt
                        %occurs in branch before), but stop pt is within
                        %current index
                    elseif start <= branchswapInd(b,1) && ...
                            stop <= branchswapInd(b,2)
                        cont_diam = [cont_diam; Branches(b). ...
                            fwhm_um(1:size(Branches(b).fwhm_um,1)- ...
                            (branchswapInd(b,2)-stop),:)];
                        %start pt is within index, as is stop pt, so use both
                    elseif start >= branchswapInd(b,1) && ...
                            stop <= branchswapInd(b,2)
                        cont_diam = [cont_diam; Branches(b). ...
                            fwhm_um(start-branchswapInd(b,1): ...
                            size(Branches(b).fwhm_um,1)- ...
                            (branchswapInd(b,2)-stop),:)];
                        %check if start and stop pts go either side of index,
                        %then use whole branch
                    elseif start <= branchswapInd(b,1) && ...
                            stop >= branchswapInd(b,2)
                        cont_diam = [cont_diam; Branches(b).fwhm_um];
                    end
                    %if branch not within range, don't do anything already
                    %created blank/filled diameter from 1st branch loop
                end %end of check if subsequent branch (2+)
            end %end of checking if branch 1
        end %end of check if first branch
        
    else
        
        %just take the diam trace from inside struct 
        cont_diam = Branches(1).fwhm_um; 
        
    end %end of checking if more than one branch to make cont trace
    
    %% save variables
    
    %output the prefs, i.e. so user can check all the parameters later
    prefs2output=prefs; 

    %save diameter, and variables needed for plotting to check code has
    %worked
    disp('saving files...'); 
    if VRflag == 0
        matfile = fullfile(exp_dir, 'contData_xyFWHM');
        save(matfile,'Branches','rawIm','calcium','stim','movement',...
            'cont_diam','time','frames','skelFlag','fps', 'prefs2output','-v7.3');
    elseif VRflag == 1
        matfile = fullfile(exp_dir, 'contData_xyFWHM');
        save(matfile,'Branches','rawIm','calcium','stim','movement',...
            'cont_diam','velocity','position','time','frames','skelFlag','fps', ...
            'prefs2output','-v7.3');
    end
    
    
end %end of looping through tif files

end %end of xyFWHM function