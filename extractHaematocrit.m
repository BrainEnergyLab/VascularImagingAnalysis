function extractHaematocrit(fname)

% created by Kira, February 2020 
% analyses RBC density (as a %) from binarised line scan data 

% input: fname
% fname is the experimental directory
% this directory should contain a tif file called 'RBCV_binary.tif' (user creates)
% this is just the velocity trace binarised - make sure it is the same num
% of frames as the loco channel - if need to remove bad frames, delete from
% both in preprocessing or after extracting 
% this directory should also contain a .ini (notepad) file, with parameter 
% info and channel 3 with the loco info

% output: continuous time series data (hct and locomotion), as well as
% extracted linescan channels (so can look for bad frames if need to
% remove)
% and saves graphs (into exp_dir) 

%% find the RBCV_binary.tif
% NB it is very important that this is binarised as it relies on that in the
% find pixels >0 
find_tif_file = findFolders(fname, 'RBCV_binary.tif');

% making approx 80ms window (as per Drew paper), must be divisble by 4
window_ms = 80; 

%for plots:
screenSz = get(0,'Screensize');

for a = 1:size(find_tif_file,2) %loop all tif files found
    
    % clear variables used throughout loop to stop interference btwn exp
    % dirs
    clearvars -except find_tif_file a fname window_ms screenSz;
    
    % disp progress to user
    disp([num2str(a),'/',num2str(size(find_tif_file,2))]);
    
    % find individual exp dir for saving vars and plots into 
    [expDir,~] = fileparts(find_tif_file{1,a});
    
    % looking for .ini file, which contains parameter info
    find_ini_file = findFolders(expDir, '*.ini');
    x = cell2mat(find_ini_file);
    ini_file = ini2struct(x);
    
    % reads the .ini file and extracts the mspline (for time vector)
    mspline = str2num(ini_file.x_.ms0x2ep0x2eline);
    fps = str2num(ini_file.x_.frames0x2ep0x2esec);
    % use mspline info to create window size
    windowsize = round((window_ms/mspline)/4)*4;

    % Get the number of frames, the width and height of the binary linescan
    info = imfinfo(find_tif_file{1,a});
    width = info(1).Width;
    height = info(1).Height;
    num_images = numel(info);
    
    % get loco channel too so can look at rest periods only
    loco_ch = cell2mat(findFolders(expDir, '*ch_3.tif'));
    % loco channel info:
    % NB num of frames and height (256) should be the same between binary
    % RBCV channel and loco channel (as this is important for timing info)
    % width can be different due to crop - so need to extract this
    % individually 
    info = imfinfo(loco_ch);
    height_ttt = info(1).Height;
    num_images_ttt = numel(info);
    width_loco = info(1).Width;
    clear info;
    
   stim_ch = cell2mat(findFolders(expDir, '*ch_4.tif'));
    
    
    %check loco and RBCV binary channel have same amount of time in them
    if height ~= height_ttt || num_images ~= num_images_ttt 
        
        disp(['Exp dir: ', expDir]);
        disp('Check num frames in RBCV vs loco ch... skipping...');
        
    else
        
        clear num_images_ttt height_ttt; 
        
        % Load in all the frames
        % inform user as process can be quite slow
        disp('loading images for all channels...');
        % Initialize an empty array to store the tiff file in.
        A = zeros(width, height, num_images);
        B = zeros(width_loco, height, num_images);
        C = zeros(width_loco, height, num_images);
        for k = 1:num_images %loop frames
            % disp every 500 frames to show progress
            if k == 1 || (mod(k,500) == 0)
                disp(['Frame ', num2str(k), '/', num2str(num_images)])
            end
            % load the binary linescan
            A(:,:,k) = imread(find_tif_file{1,a}, k)';
            % load locomotion
            B(:,:,k) = imread(loco_ch, k)';
            %load stim
            C(:, :, k) = imread(stim_ch, k)';
        end %end of frames loop
        
        % reorder the data:
        % put all frames on top of each other, so it becomes a 2D image
        % get the time x space line scan info
        % binary image:
        binaryLine = reshape(A, [width, height*num_images]);
        locoLine = reshape(B, [width_loco, height*num_images]);
        stimLine = reshape(C, [width_loco, height*num_images]);
        %dnt care about width of loco channel, just take avg for walking or not
        locoLine = nanmean(locoLine,1);
        stimLine = nanmean(stimLine, 1);
        clear A B C;
        
        %% locomotion - extract as time series - avg loco per time pt
        
        % call function:
        % get out avg loco from ch3 for each point of the sliding time window
        [~,locomotion] = avgByTimeWindow(locoLine',windowsize,0);
        [~,stim] = avgByTimeWindow(stimLine',windowsize,0);
        % clean loco so below thresh removed, and data is between 0 and 1
        % so can detect loco events off and on later
        [locomotion] = cleanLoco(locomotion);
        [locomotion] = norm_01(locomotion);
        locomotion(:,find(locomotion<=0.01))=0;
        
         %find when stim is on, binarise stim channel:
        stim = stim>(mean(stim)+std(stim));
        
        %call function to check if any spikes in stim signal (which last
        %less than 1s)
        %it will replace the spikes with zeros - this must be a binarised
        %stim trace
        [stim] = findRogueStims(stim',fps);
        
        %% haematocrit
        
        % call function
        % calculate haematocrit % for multiple time points along data
        [time_ttt,hct] = avgByTimeWindow(binaryLine',windowsize,1);
        if size(hct,1) > size(hct,2)
            hct = hct';
        end
        
        % create time vector for plots
        % convert to be in seconds (using mspline info)
        time = (time_ttt*mspline)/1000;
        if size(time,1) > size(time,2)
            time = time';
        end
        clear time_ttt;
        
        %% plots:
        
        % plot the original line scan continuous data, with thresholded data
        % i.e. check the threshold has worked - as this will be used for % calc
        figure;
        %make fig size of screen
        set(gcf, 'Position', [screenSz(1) screenSz(2) screenSz(3) ...
            screenSz(4)]);
        a1=subplot(3,1,1);
        imagesc(time,[],binaryLine);
        colormap gray;
        title('Original Binary Line Scan Data');
        xlabel('Time (s)');
        a2=subplot(3,1,2);
        plot(time, hct);
        title('Haematocrit Plot');
        xlabel('Time (s)');
        ylabel('% RBC density');
        a3=subplot(3,1,3);
        plot(time, locomotion);
        title('Locomotion Trace');
        xlabel('Time (s)');
        ylabel('A.U.');
        linkaxes([a1,a2,a3],'x');
        figSave = 'HctTrace.png';
        saveas(gcf, fullfile([expDir,filesep,figSave]));
        close;
        
        %% save vars
        save([expDir,filesep,'contData_ls_Hct'],'time','hct',...
            'binaryLine', 'locoLine', 'locomotion', 'stim');
        
    end
    
end %end of loop tif files

end %end of func

