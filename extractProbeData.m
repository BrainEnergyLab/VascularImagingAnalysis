function extractProbeData(fname)

%function written by Kira, August 2017, updated May 2018
%updated Nov 2020
%function searches through experimental directories for excel files
%containing oxyprobe/LDF data 
%it will then extract the data as a .mat file, and save into the directory
%
%
%INPUTS-
%fname = directory which contaisn excel sheets 
%NB/ currently this will search for ALL excel sheets in the folders
%inputted, to not limit name choice, so no other excel sheets should be in
%the folders inputted or it will error 
%can update and force user to have 'spont' or 'stim' in the name if needed
%OUTPUTS-
%no variables are outputted, a .mat file, with the data from the excel now
%as a matlab file will be saved into the local exp dirs
%
%other functions needed in path to run this function:
%findFolders

%search for all excel sheets within the specified directory 
find_excel_file = findFolders(fname, '*.xlsx');

%loop to open each excel file
for i = 1:size(find_excel_file,2)
    
    %inform user of progress
    disp(['processing file ',num2str(i),'/',num2str(size( ...
        find_excel_file,2))]);
    
    %read the probe data in as a table with headings (in case the variables
    %switch places)
    data_ttt = readtable(find_excel_file{1,i});
    %create a blank matrix to fill the data from the table into in the
    %correct order
    %1=flux, 2=so2, 3=speed, 4=hbo, 5=hbr, 6=hbt, 7=cmro2
    data = zeros(7,size(data_ttt,1)); 
    
    %get column names - safer to put data into right place using name
    colNames = data_ttt.Properties.VariableNames;
    %loop column names, and put haem data into associated column
    % i.e. %1=flux, 2=so2, 3=speed, 4=hbo, 5=hbr, 6=hbt, 7=cmro2
    for a = 1:size(colNames,2)
        title = colNames{1,a};
        if strmatch(title, 'Flux', 'exact')
            data(1,:) = table2array(data_ttt(:,a));
        elseif strmatch(title, 'SO2', 'exact')
            data(2,:) = table2array(data_ttt(:,a));
        elseif strmatch(title, 'Speed', 'exact')
            data(3,:) = table2array(data_ttt(:,a));
        elseif strmatch(title, 'oxyHb', 'exact')
            data(4,:) = table2array(data_ttt(:,a));
        elseif strmatch(title, 'deoxyHb', 'exact')
            data(5,:) = table2array(data_ttt(:,a));
        elseif strmatch(title, 'totalHb', 'exact')
            data(6,:) = table2array(data_ttt(:,a));
        end
        clear title;
    end

    %calculate CMRO2
    %(flux*hbr)/hbt
    data(7,:) = (data(1,:)).*(data(5,:)./data(6,:));
    
    %loco info is in Analogue 1
    %this trace will look different depending if readRotEncoder or virmen 
    %was used during recording (this will be sorted in later code) -- 
    %for now just extract channel
    Index = find(contains(colNames,'Ana1'));
    movement = table2array(data_ttt(:,Index))';
    clear Index; 
    %process raw loco data
    [locomotion] = cleanLoco(movement);
    locomotion_raw = abs(diff(movement));
    locomotion_raw = [locomotion_raw,0];
    
    %stim info is in Analogue 2
    Index = find(contains(colNames,'Ana2'));
    stim = table2array(data_ttt(:,Index))';
    clear Index; 
    %find when stim is on, binarise stim channel:
    stim = stim > (mean(stim)+std(stim));
    %set first 50 points to zero, in case of initial spike
    stim(1:50) = 0;
    
    %create frame vector (for plots)
    frames = [1:size(data,2)];
    
    %fps of Hb-LDF probe is 40
    fps = 40; %40Hz
    %create time vector (for plots)
    time = frames/fps;
    
    %save data into correct exp_dir
    [exp_dir,~] = fileparts(find_excel_file{1,i});
    matfile = fullfile(exp_dir, 'contData');
    save(matfile,'data','movement','locomotion','locomotion_raw',...
        'stim','time','frames','fps');
    
end %end of loop excel files

end %end of function
