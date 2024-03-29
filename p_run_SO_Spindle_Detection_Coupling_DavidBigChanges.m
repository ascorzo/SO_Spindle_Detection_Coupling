%% Define some important variables
s_TotalTimeSec = 30; %total time in seconds we are analysing
s_fs = 1000; % sample frequency
s_TimeSam = s_TotalTimeSec * s_fs; %total time in samples

% In case the time is unbalanced before and after the stimuli, specify the times.
s_TimeAfterCero = 15;
s_TimeBeforeCero = 15;

%% Set up user land
% Normally not needed to add slash to end of folder paths. But necessary in
% case of searching for specific extensions.
if contains(computer,'PCWIN')
    slashSys = '\';
else
    slashSys = '/';
end

%Selection of folder where channel data is contained
ChanVsSource = questdlg('Import channel recording data or source-computed cortex data?', ...
    'Channel or Source?', ...
    'Channel','Source','Channel');

%Selection of folder where source data is contained
% OnVsOff = questdlg('Select if you want to evaluate off condition or odor condition', ...
%     'Odor On or Odor Off?', ...
%     'On','Off','On');


%% Create a Files List in order to go through them for the analysis

switch ChanVsSource
    case 'Channel'
        %Containing channel data
        if ~exist('pathNameChan','var') || size(pathNameChan,2) < 3
            % path to datasets containing channel data
            pathNameChan = [uigetdir(cd,'Choose the folder that contains the CHANNEL datasets'), slashSys];
            addpath(pathNameChan)
            FilesListChanOdor = dir([pathNameChan,'*Odor.mat']);
            FilesListChanPlacebo = dir([pathNameChan,'*Placebo.mat']);
            
            FilesList = FilesListChanOdor;
        end
    case 'Source'
        if  ~exist('pathNameSource','var') || size(pathNameSource,2) < 3
            %path to datasets with source estimation
            pathNameSource = [uigetdir(cd,'Choose the folder that contains the SOURCE datasets'), slashSys];
            addpath(pathNameSource)
            FilesListSourceOdor = dir([pathNameSource,'*Odor.mat']);
            FilesListSourcePlacebo = dir([pathNameSource,'*Placebo.mat']);
            
            FilesList = FilesListSourceOdor;
        end
end

%% Select Which areas to compare

if strcmp(ChanVsSource,'Channel')
    FileTemp = load(FilesListChanOdor(1).name);
    Sources = FileTemp.Channel.label;
    %Asks for selection of a channel to detect the Slow Oscillations and
    %assigns to str_ChanSO
    [indx,tf] = listdlg('PromptString','Select a channel: for Slow Oscillations',...
        'SelectionMode','single','ListSize',[400 400],...
        'ListString',Sources);
    str_ChanSO = Sources(indx);
    
    %Asks for selection of a channel to detect the Sleep Spindles and
    %assigns to str_ChanSS
    [indx,tf] = listdlg('PromptString','Select a channel: for Sleep Spindles',...
        'SelectionMode','single','ListSize',[400 400],...
        'ListString',Sources);
    str_ChanSS = Sources(indx);
    
else
    FileTemp = load(FilesListSourceOdor(1).name);
    Sources = string({FileTemp.Atlas.Scouts.Label});
    %Asks for selection of a region to detect the Slow Oscillations and
    %assigns to str_ROI_SO
    [indx,tf] = listdlg('PromptString','Select a source for Slow Oscillations:',...
        'SelectionMode','single','ListSize',[400 400],...
        'ListString',Sources);
    str_ROI_SO = Sources(indx);
    %Asks for selection of a region to detect the Sleep Spindles and
    %assigns to str_ROI_SS
    [indx,tf] = listdlg('PromptString','Select a source for Sleep Spindles:',...
        'SelectionMode','single','ListSize',[400 400],...
        'ListString',Sources);
    str_ROI_SS = Sources(indx);
end

%% Loading Odor sets into memory

% Speed up script: Load files only when needed after changing user input of
% type of data to compute
if ~exist('subjectFilesOdor','var') ||...
        ~numel(FilesList) == numel(subjectFilesOdor) ||...
        ~exist('previous_ChanVsSource','var') ||...
        ~strcmp(previous_ChanVsSource, ChanVsSource)
    for Load2Mem = 1:numel(FilesList)
        if strcmp(ChanVsSource,'Channel')
            subjectFilesOdor{Load2Mem,1} = load([pathNameChan FilesListChanOdor(Load2Mem).name]);
        else
            % subjectFilesOdor{Load2Mem,1} = load([pathNameSource FilesListSourceOdor(Load2Mem).name]);
            matfiles_subjectFilesOdor{Load2Mem,1} = matfile([pathNameSource FilesListSourceOdor(Load2Mem).name]);
            
            % Loads everything except 'Value' field. This field is the
            % heavy part of the files. Later, only needed data will be
            % loeded.
            subjectFilesOdor{Load2Mem,1} = load(matfiles_subjectFilesOdor{Load2Mem,1}.Properties.Source,...
                '-regexp', '^(?!Value)\w');
        end
    end
end

% Append "Value" field to subject files with only the data of interest
for Load2Mem = 1:numel(FilesList)
    if strcmp(ChanVsSource,'Source')
        
        % Determine line of 'Value' to load from subject file in loop below
        s_ROI_SO = strcmp({subjectFilesOdor{Load2Mem, 1}.Atlas.Scouts.Label},str_ROI_SO);
        idx_ROI_SO = find(s_ROI_SO);
        s_ROI_SS = strcmp({subjectFilesOdor{Load2Mem, 1}.Atlas.Scouts.Label},str_ROI_SS);
        idx_ROI_SS = find(s_ROI_SS);
        
        % line 1 will be source of SO, line 2 source of SS
        subjectFilesOdor{Load2Mem,1}.Value(1,:) = matfiles_subjectFilesOdor{Load2Mem,1}.Value(idx_ROI_SO,:);
        subjectFilesOdor{Load2Mem,1}.Value(2,:) = matfiles_subjectFilesOdor{Load2Mem,1}.Value(idx_ROI_SS,:);
        
    end
end

%% Loading Placebo sets into memory

if ~exist('subjectFilesPlacebo','var') ||...
        ~numel(FilesList) == numel(subjectFilesPlacebo) ||...
        ~exist('previous_ChanVsSource','var') ||...
        ~strcmp(previous_ChanVsSource, ChanVsSource)
    for Load2Mem = 1:numel(FilesList)
        if strcmp(ChanVsSource,'Channel')
            subjectFilesPlacebo{Load2Mem,1} = load([pathNameChan FilesListChanPlacebo(Load2Mem).name]);
        else
            matfiles_subjectFilesPlacebo{Load2Mem,1} = matfile([pathNameSource FilesListSourcePlacebo(Load2Mem).name]);
            
            % Again, loads everything except 'Value' field
            subjectFilesPlacebo{Load2Mem,1} = load(matfiles_subjectFilesPlacebo{Load2Mem,1}.Properties.Source,...
                '-regexp', '^(?!Value)\w');
        end
    end
end

% Append "Value" field to subject files with only the data of interest
for Load2Mem = 1:numel(FilesList)
    if strcmp(ChanVsSource,'Source')
        
        % Determine line of 'Value' to load from subject file in loop below
        s_ROI_SO = strcmp({subjectFilesPlacebo{Load2Mem, 1}.Atlas.Scouts.Label},str_ROI_SO);
        idx_ROI_SO = find(s_ROI_SO);
        s_ROI_SS = strcmp({subjectFilesPlacebo{Load2Mem, 1}.Atlas.Scouts.Label},str_ROI_SS);
        idx_ROI_SS = find(s_ROI_SS);
        
        % line 1 will be source of SO, line 2 source of SS
        subjectFilesPlacebo{Load2Mem,1}.Value(1,:) = matfiles_subjectFilesPlacebo{Load2Mem,1}.Value(idx_ROI_SO,:);
        subjectFilesPlacebo{Load2Mem,1}.Value(2,:) = matfiles_subjectFilesPlacebo{Load2Mem,1}.Value(idx_ROI_SS,:);
    end
end

%% Going first through Off, then through On period
for run = 1:2
    if run == 1
        OnVsOff = 'Off';
    elseif run == 2
        OnVsOff = 'On';
    end
    
    fprintf('<!> Going now through "%s" period\n',OnVsOff);
    
    %% Odor datasets
    
    for subj = 1:length(subjectFilesOdor)
        
        switch ChanVsSource
            case 'Channel'
                %Set up datasets containing channel data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fileChan = load(FilesListChanOdor(subj).name);
                s_NumChannels = length(fileChan.Channel.number);
                v_Time = fileChan.Channel.times;
                s_Trials = fileChan.Channel.trials;
                DataOdorChan = double(fileChan.Channel.data);
                
            case 'Source'
                %Set up datasets containing source data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fileSource = subjectFilesOdor{subj};
                s_NumScouts = size(fileSource.Value,1); % Detect the number of scouts
                v_Time = fileSource.Time(1:s_TimeSam);
                s_Trials = size(fileSource.Value,2)/s_TimeSam;
                DataOdorSource = reshape(fileSource.Value,[s_NumScouts,s_TimeSam,s_Trials]);
        end
        
        
        for trial = 1:s_Trials
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  SO detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %-----------------------------------------
            % Select Data from where to detect the SO
            %-----------------------------------------
            switch ChanVsSource
                case 'Channel'
                    findChannel = find(strcmp(fileChan.Channel.label, str_ChanSO));
                    DataOdorTrial_SO = DataOdorChan(findChannel,:,trial);
                case 'Source'
                    % s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                    DataOdorTrial_SO = DataOdorSource(1,:,trial); % 1,y,z because s_ROI_SO has been put in first line
            end
            
            %------------------------------------------------
            % Definition of variables to do the SO detection
            %------------------------------------------------
            Label = 'SO';
            SleepScoring = ones(length(DataOdorTrial_SO),1);
            ScoringArtefacts = zeros(length(DataOdorTrial_SO),1);
            if strcmp(OnVsOff,'On') == 1
                % For the 'On' Period, we assign the values before the trigger
                % is presented as artefacts.
                ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
            else
                % For the 'Off' Period, we assign the values after the trigger
                % is presented as artefacts.
                ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
            end
            
            %---------------
            % SO detection
            %---------------
            [FilteredEEG_SO, out_SO] = SODetection(DataOdorTrial_SO', SleepScoring, ScoringArtefacts, Label, s_fs);
            %FilteredEEG_SO contains the filtered signal, and out_SO contains
            %all the information related to the SO detection
            
            v_time = (0:1/s_fs:(length(FilteredEEG_SO)-1)/s_fs) - s_TimeBeforeCero;
            
            %----------------------------------------------
            % Plot in case you want to check the detection
            %----------------------------------------------
            subplot(2,1,1);
            plot(v_time,FilteredEEG_SO);
            for i=1:size(out_SO.trialinfo,1)
                line([out_SO.trialinfo(i,2)/out_SO.fsample out_SO.trialinfo(i,4)/out_SO.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
            end
            %ylim([-.3 .3])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Spindle detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Store choosen region for later saving
            LabelCortex = Label;
            
            %-----------------------------------------
            % Select Data from where to detect the SS
            %-----------------------------------------
            switch ChanVsSource
                case 'Channel'
                    findChannel = find(strcmp(fileChan.Channel.label, str_ChanSS));
                    DataOdorTrial_SS = DataOdorChan(findChannel,:,trial);
                case 'Source'
                    % s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                    DataOdorTrial_SS = DataOdorSource(2,:,trial); % 2,y,z because s_ROI_SS has been put in first line
            end
            
            %-------------------------------------------------
            % Definition of variables to do the SS detection
            %-------------------------------------------------
            Label = 'SS';
            SleepScoring = ones(length(DataOdorTrial_SS),1);
            ScoringArtefacts = zeros(length(DataOdorTrial_SS),1);
            
            if strcmp(OnVsOff,'On') == 1
                % For the 'On' Period, we assign the values before the trigger
                % is presented as artefacts.
                ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
            else
                % For the 'Off' Period, we assign the values after the trigger
                % is presented as artefacts.
                ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
            end
            
            %---------------
            % SS detection
            %---------------
            [FilteredEEG_SS, out_SS] = SpindleDetection(DataOdorTrial_SS', SleepScoring, ScoringArtefacts, Label, s_fs);
            %FilteredEEG_SS contains the filtered signal, and out_SS contains
            %all the information related to the SS detection
            
            %----------------------------------------------
            % Plot in case you want to check the detection
            %----------------------------------------------
            subplot(2,1,2);
            plot(v_time,FilteredEEG_SS);
            for i=1:size(out_SS.trialinfo,1)
                line([out_SS.trialinfo(i,2)/out_SS.fsample out_SS.trialinfo(i,4)/out_SS.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
            end
            %ylim([-.3 .3])
            
            %-------------------------------------------------------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Phase Coupling calculation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Calculate the Density of spindles in each Bin of SO phase
            mode = 'Density';
            [spindleInBininTrial(trial,:),xedges_Odor] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
            
            %Calculate the Power of spindles in each Bin of SO phase
            mode = 'Power';
            [spindlePowerInBininTrial(trial,:),xedges_Odor] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
            
            % Calculates the total density of spindles
            v_DensitySpindles(trial) = size(out_SS.trialinfo,1);
            % Calculates the mean of the power of all spindles
            v_MeanPowerSpindles(trial) = mean(out_SS.trialinfo(:,12));
            
            %Create Total Vectors in order to save variables at the end
            spindleInBinTotalOdor.(strcat('S',num2str(subj))).(strcat('T',num2str(trial)))= spindleInBininTrial(trial,:);
            spindlePowerInBinTotalOdor.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = spindlePowerInBininTrial(trial,:);
            DensitySpindlesTotalOdor.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = v_DensitySpindles(trial);
            PowerSpindlesTotalOdor.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = out_SS.trialinfo(:,12);
        end
        
        v_MeanDensitySpindlesOdor(subj) = nanmean(v_DensitySpindles);
        v_MeanPowerSpindlesOdor(subj) = nanmean(v_MeanPowerSpindles);
        spindleInBininSubject_Odor(subj,:) = nanmean(spindleInBininTrial,1);
        spindlePowerInBininSubject_Odor(subj,:) = nanmean(spindlePowerInBininTrial,1);
    end
    
    %% Placebo datasets
    for subj = 1:length(subjectFilesPlacebo)
        
        switch ChanVsSource
            case 'Channel'
                %Set up datasets containing channel data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fileChan = load(FilesListChanPlacebo(subj).name);
                s_NumChannels = length(fileChan.Channel.number);
                v_Time = fileChan.Channel.times;
                s_Trials = fileChan.Channel.trials;
                DataPlaceboChan = double(fileChan.Channel.data);
                
            case 'Source'
                %Set up datasets containing source data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fileSource = subjectFilesPlacebo{subj};
                s_NumScouts = size(fileSource.Value,1); % Detect the number of scouts
                v_Time = fileSource.Time(1:s_TimeSam);
                s_Trials = size(fileSource.Value,2)/s_TimeSam;
                DataPlaceboSource = reshape(fileSource.Value,[s_NumScouts,s_TimeSam,s_Trials]);
                
        end
        
        for trial = 1:s_Trials
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  SO detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %-----------------------------------------
            % Select Data from where to detect the SO
            %-----------------------------------------
            switch ChanVsSource
                case 'Channel'
                    findChannel = find(strcmp(fileChan.Channel.label, str_ChanSO));
                    DataPlaceboTrial_SO = DataPlaceboChan(findChannel,:,trial);
                case 'Source'
                    % s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                    DataPlaceboTrial_SO = DataPlaceboSource(1,:,trial);
            end
            
            %------------------------------------------------
            % Definition of variables to do the SO detection
            %------------------------------------------------
            Label = 'SO';
            SleepScoring = ones(length(DataPlaceboTrial_SO),1);
            ScoringArtefacts = zeros(length(DataPlaceboTrial_SO),1);
            if strcmp(OnVsOff,'On') == 1
                % For the 'On' Period, we assign the values before the trigger
                % is presented as artefacts.
                ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
            else
                % For the 'Off' Period, we assign the values after the trigger
                % is presented as artefacts.
                ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
            end
            
            [FilteredEEG_SO, out_SO] = SODetection(DataPlaceboTrial_SO', SleepScoring, ScoringArtefacts, Label, s_fs);
            %FilteredEEG_SO contains the filtered signal, and out_SO contains
            %all the information related to the SO detection
            
            v_time = (0:1/s_fs:(length(FilteredEEG_SO)-1)/s_fs) - s_TimeBeforeCero;
            
            %------------------------------------------------------------------
            % Plot in case you want to check the detection
            %------------------------------------------------------------------
            subplot(2,1,1);
            plot(v_time,FilteredEEG_SO);
            for i=1:size(out_SO.trialinfo,1)
                line([out_SO.trialinfo(i,2)/out_SO.fsample out_SO.trialinfo(i,4)/out_SO.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
            end
            %ylim([-.3 .3])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Spindle detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Store choosen region for later saving
            LabelCortex = Label;
            
            %-----------------------------------------
            % Select Data from where to detect the SS
            %-----------------------------------------
            switch ChanVsSource
                case 'Channel'
                    findChannel = find(strcmp(fileChan.Channel.label, str_ChanSS));
                    DataPlaceboTrial_SS = DataPlaceboChan(findChannel,:,trial);
                case 'Source'
                    % s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                    DataPlaceboTrial_SS = DataPlaceboSource(2,:,trial);
            end
            
            %-------------------------------------------------
            % Definition of variables to do the SS detection
            %-------------------------------------------------
            Label = 'SS';
            SleepScoring = ones(length(DataPlaceboTrial_SS),1);
            ScoringArtefacts = zeros(length(DataPlaceboTrial_SS),1);
            if strcmp(OnVsOff,'On') == 1
                % For the 'On' Period, we assign the values before the trigger
                % is presented as artefacts.
                ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
            else
                % For the 'Off' Period, we assign the values after the trigger
                % is presented as artefacts.
                ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
            end
            
            [FilteredEEG_SS, out_SS] = SpindleDetection(DataPlaceboTrial_SS', SleepScoring, ScoringArtefacts, Label, s_fs);
            %FilteredEEG_SS contains the filtered signal, and out_SS contains
            %all the information related to the SS detection
            
            %----------------------------------------------
            % Plot in case you want to check the detection
            %----------------------------------------------
            subplot(2,1,2);
            plot(v_time,FilteredEEG_SS);
            for i=1:size(out_SS.trialinfo,1)
                line([out_SS.trialinfo(i,2)/out_SS.fsample out_SS.trialinfo(i,4)/out_SS.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
            end
            %ylim([-.3 .3])
            %-------------------------------------------------------------------------
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Phase Coupling calculation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Calculate the Density of spindles in each Bin of SO phase
            mode = 'Density';
            [spindleInBininTrial(trial,:),xedges_Placebo] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
            
            %Calculate the Power of spindles in each Bin of SO phase
            mode = 'Power';
            [spindlePowerInBininTrial(trial,:),~] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
            
            % Calculates the total density of spindles
            v_DensitySpindles(trial) = size(out_SS.trialinfo,1);%/s_TimeAfterCero; % Play with density or just number (Density as # of spindles in time)
            % Calculates the mean of the power of all spindles
            v_MeanPowerSpindles(trial) = mean(out_SS.trialinfo(:,12));
            
            %Create Total Vectors in order to save variables at the end
            spindleInBinTotalPlacebo.(strcat('S',num2str(subj))).(strcat('T',num2str(trial)))= spindleInBininTrial(trial,:);
            spindlePowerInBinTotalPlacebo.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = spindlePowerInBininTrial(trial,:);
            DensitySpindlesTotalPlacebo.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = v_DensitySpindles(trial);
            PowerSpindlesTotalPlacebo.(strcat('S',num2str(subj))).(strcat('T',num2str(trial))) = out_SS.trialinfo(:,12);
        end
        v_MeanDensitySpindlesPlacebo(subj) = nanmean(v_DensitySpindles);
        v_MeanPowerSpindlesPlacebo(subj) = nanmean(v_MeanPowerSpindles);
        spindleInBininSubject_Placebo(subj,:) = nanmean(spindleInBininTrial,1);
        spindlePowerInBininSubject_Placebo(subj,:) = nanmean(spindlePowerInBininTrial,1);
    end
    
    %% Define variables for figure titles later
    if strcmp(ChanVsSource,'Channel')
        SOsource = str_ChanSO;
        Spindlesource = str_ChanSS;
    else
        SOsource = str_ROI_SO;
        Spindlesource = str_ROI_SS;
    end
    
    %% Plots
    % figure
    % histogram(v_MeanDensitySpindlesOdor)
    % hold on
    % histogram(v_MeanDensitySpindlesPlacebo)
    % title('Density')
    
    
    % figure
    % histogram(v_MeanPowerSpindlesOdor)
    % hold on
    % histogram(v_MeanPowerSpindlesPlacebo)
    % title('Power')
    
    
    %% Polar Plot for Density
    meanSpindleInBin_Placebo = nanmean(spindleInBininSubject_Placebo,1);
    meanSpindleInBin_Odor = nanmean(spindleInBininSubject_Odor,1);
    
    spindleInBininSubject.Placebo = spindleInBininSubject_Placebo;
    spindleInBininSubject.Odor = spindleInBininSubject_Odor;
    save([cd, 'MeanSpindlesInBin'],'spindleInBininSubject');
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    figNumber = polarhistogram('BinEdges', xedges_Placebo, 'BinCounts', meanSpindleInBin_Placebo,'FaceColor', 'blue',...
        'FaceAlpha',.2);
    hold on;
    figNumber = polarhistogram('BinEdges', xedges_Odor, 'BinCounts', meanSpindleInBin_Odor,'FaceColor', 'red',...
        'FaceAlpha',.2);
    thetaticks(0:90:315); %rlim([0 1.2]);
    pax = gca; pax.ThetaAxisUnits = 'radians';
    pax.FontSize = 12; pax.GridColor = 'black';
    legend('sham','cue');
    title(strcat('Spindle density of', {' '}, Spindlesource,...
        ' phase-locked to Slow Osc. of', {' '}, SOsource));
    
    %% Polar Plot for Power
    
    meanSpindlePowerInBin_Placebo = nanmean(spindlePowerInBininSubject_Placebo,1);
    meanSpindlePowerInBin_Odor = nanmean(spindlePowerInBininSubject_Odor,1);
    
    spindlePowerInBininSubject.Placebo = spindlePowerInBininSubject_Placebo;
    spindlePowerInBininSubject.Odor = spindlePowerInBininSubject_Odor;
    save([cd, 'MeanSpindlesPowerInBin'],'spindlePowerInBininSubject');
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    figPower = polarhistogram('BinEdges', xedges_Placebo, 'BinCounts', meanSpindlePowerInBin_Placebo,'FaceColor', 'blue',...
        'FaceAlpha',.2);
    hold on;
    figPower = polarhistogram('BinEdges', xedges_Odor, 'BinCounts', meanSpindlePowerInBin_Odor,'FaceColor', 'red',...
        'FaceAlpha',.2);
    thetaticks(0:90:315); %rlim([0 2.2]);
    pax = gca; pax.ThetaAxisUnits = 'radians';
    pax.FontSize = 12; pax.GridColor = 'black';
    legend('sham','cue');
    title(strcat('Spindle power of', {' '}, Spindlesource,...
        ' phase-locked to Slow Osc. of', {' '}, SOsource));
    
    %% Save useful information
    % =========================================================================
    % Be aware that changing the way matrices holding results are generated
    % will greatly affect posterior scripts of table generation and statistical
    % computation and will most likely create breakage of the scripts!
    % Please only do changes if really needed.
    % =========================================================================
    
    % Prepare path for saving data
    savePath = strcat(cd, slashSys, 'Results', slashSys);
    
    % Prepare structure that holds all results
    if strcmp(ChanVsSource,'Channel')
        Channel.LabelSO = str_ChanSO;
        Channel.LabelSpindle = str_ChanSS;
        Channel.OnVsOff = OnVsOff;
        Channel.spindleInBinTotalOdor = spindleInBinTotalOdor;
        Channel.spindleInBinTotalPlacebo = spindleInBinTotalPlacebo;
        Channel.spindlePowerInBinTotalOdor = spindlePowerInBinTotalOdor;
        Channel.spindlePowerInBinTotalPlacebo = spindlePowerInBinTotalPlacebo;
        Channel.DensitySpindlesTotalOdor = DensitySpindlesTotalOdor;
        Channel.DensitySpindlesTotalPlacebo = DensitySpindlesTotalPlacebo;
        Channel.PowerSpindlesTotalOdor = PowerSpindlesTotalOdor;
        Channel.PowerSpindlesTotalPlacebo = PowerSpindlesTotalPlacebo;
        
    elseif strcmp(ChanVsSource,'Source')
        Source.LabelSO = str_ROI_SO;
        Source.LabelSpindle = str_ROI_SS;
        Source.OnVsOff = OnVsOff;
        Source.spindleInBinTotalOdor = spindleInBinTotalOdor;
        Source.spindleInBinTotalPlacebo = spindleInBinTotalPlacebo;
        Source.spindlePowerInBinTotalOdor = spindlePowerInBinTotalOdor;
        Source.spindlePowerInBinTotalPlacebo = spindlePowerInBinTotalPlacebo;
        Source.DensitySpindlesTotalOdor = DensitySpindlesTotalOdor;
        Source.DensitySpindlesTotalPlacebo = DensitySpindlesTotalPlacebo;
        Source.PowerSpindlesTotalOdor = PowerSpindlesTotalOdor;
        Source.PowerSpindlesTotalPlacebo = PowerSpindlesTotalPlacebo;
        
    end
    
    % Info for plotting polar histogram
    PlotInfo.meanSpindleInBin_Odor = meanSpindleInBin_Odor;
    PlotInfo.meanSpindleInBin_Placebo = meanSpindleInBin_Placebo;
    PlotInfo.meanSpindlePowerInBin_Odor = meanSpindlePowerInBin_Odor;
    PlotInfo.meanSpindlePowerInBin_Placebo = meanSpindlePowerInBin_Placebo;
    PlotInfo.xedges_Odor = xedges_Odor;
    PlotInfo.xedges_Placebo = xedges_Placebo;
    
    % Save results
    if ~exist(strcat(cd, slashSys, 'Results'),'dir')
        mkdir (strcat(cd, slashSys, 'Results'));
    end
    
    if strcmp(ChanVsSource,'Channel')
        save(strcat(savePath,str_ChanSO,'_vs_',str_ChanSS,'_',OnVsOff,'.mat'), 'Channel', 'PlotInfo');
        % Save figure
        saveas(figNumber, strcat(savePath,str_ChanSO,'_vs_',str_ChanSS,'_',OnVsOff,'_Density.bmp'));
        saveas(figPower, strcat(savePath,str_ChanSO,'_vs_',str_ChanSS,'_',OnVsOff,'_Power.bmp'));
    else
        save(strcat(savePath,str_ROI_SO,'_vs_',str_ROI_SS,'_',OnVsOff,'.mat'), 'Source', 'PlotInfo');
        % Save figure
        saveas(figNumber, strcat(savePath,str_ROI_SO,'_vs_',str_ROI_SS,'_',OnVsOff,'_Density.bmp'));
        saveas(figPower, strcat(savePath,str_ROI_SO,'_vs_',str_ROI_SS,'_',OnVsOff,'_Power.bmp'));
    end
    
end

% Important in order to speed up script
previous_ChanVsSource = ChanVsSource;
