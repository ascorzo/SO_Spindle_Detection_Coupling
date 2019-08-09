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
OnVsOff = questdlg('Select if you want to evaluate off condition or odor condition', ...
    'Odor On or Odor Off?', ...
    'On','Off','On');


%% Create a Files List in order to go through them for the analysis

switch ChanVsSource
    case 'Channel'
        %Containing channel data
        if ~exist('pathNameChan','var')
            % path to datasets containing channel data
            pathNameChan = [uigetdir(cd,'Choose the folder that contains the CHANNEL datasets'), slashSys];
            addpath(pathNameChan)
            FilesListChanOdor = dir([pathNameChan,'*Odor.mat']);
            FilesListChanPlacebo = dir([pathNameChan,'*Placebo.mat']);
            
            FilesList = FilesListChanOdor;
        end
    case 'Source'
        if ~exist('pathNameSource','var')
            %path to datasets with source estimation
            pathNameSource = [uigetdir(cd,'Choose the folder that contains the SOURCE datasets'), slashSys];
            addpath(pathNameSource)
            FilesListSourceOdor = dir([pathNameSource,'*Odor.mat']);
            FilesListSourcePlacebo = dir([pathNameSource,'*Placebo.mat']);
            
            FilesList = FilesListSourceOdor;
        end
end

%% Loading Odor sets into memory

for Load2Mem = 1:numel(FilesList)
    if strcmp(ChanVsSource,'Channel')
        subjectFiles{Load2Mem,1} = load([pathNameChan FilesListChanOdor(Load2Mem).name]);
    else
        subjectFiles{Load2Mem,1} = load([pathNameSource FilesListSourceOdor(Load2Mem).name]);
    end
end
 

%% Select Which areas to compare 

if strcmp(ChanVsSource,'Channel')
    FileTemp = load(FilesListChanOdor(1).name);
    Sources = FileTemp.Channel.label;
    %Asks for selection of a channel to detect the Slow Oscillations and
    %assigns to str_ChanSO
    [indx,tf] = listdlg('PromptString','Select a source: for Slow Oscillations',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ChanSO = Sources(indx);
    
    %Asks for selection of a channel to detect the Sleep Spindles and
    %assigns to str_ChanSS
    [indx,tf] = listdlg('PromptString','Select a source: for Sleep Spindles',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ChanSS = Sources(indx);
    
else
    FileTemp = subjectFiles{1};
    Sources = string({FileTemp.Atlas.Scouts.Label});
    %Asks for selection of a region to detect the Slow Oscillations and
    %assigns to str_ROI_SO
    [indx,tf] = listdlg('PromptString','Select a source for Slow Oscillations:',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ROI_SO = Sources(indx);
    %Asks for selection of a region to detect the Sleep Spindles and
    %assigns to str_ROI_SS
    [indx,tf] = listdlg('PromptString','Select a source for Sleep Spindles:',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ROI_SS = Sources(indx);
end

%% Odor datasets

for subj = 1:length(subjectFiles)
    
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
            fileSource = subjectFiles{subj};
            s_NumScouts = length(fileSource.Atlas.Scouts); % Detect the number of scouts
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
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                DataOdorTrial_SO = DataOdorSource(s_ROI,:,trial);
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
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                DataOdorTrial_SS = DataOdorSource(s_ROI,:,trial);
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
        [spindleInBininTrial(trial,:),xedges_dor] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
        
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
    
    v_MeanDensitySpindlesOdor(subj) = nanmean(v_DensitySpindles); %you can do either mean or sum
    v_MeanPowerSpindlesOdor(subj) = nanmean(v_MeanPowerSpindles); %you can do either mean or sum
    spindleInBininSubject_Odor(subj,:) = nanmean(spindleInBininTrial,1); %you can do either mean or sum
    spindlePowerInBininSubject_Odor(subj,:) = nanmean(spindlePowerInBininTrial,1); %you can do either mean or sum
end

clear subjectFilesSource
%% Loading Placebo sets into memory

for Load2Mem = 1:numel(FilesList)
    if strcmp(ChanVsSource,'Channel')
        subjectFiles{Load2Mem,1} = load([pathNameChan FilesListChanPlacebo(Load2Mem).name]);
    else
        subjectFiles{Load2Mem,1} = load([pathNameSource FilesListSourcePlacebo(Load2Mem).name]);
    end
end

%% Placebo datasets
for subj = 1:length(subjectFiles)
    
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
            fileSource = subjectFiles{subj};
            s_NumScouts = length(fileSource.Atlas.Scouts); % Detect the number of scouts
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
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                DataPlaceboTrial_SO = DataPlaceboSource(s_ROI,:,trial);
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
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                DataPlaceboTrial_SS = DataPlaceboSource(s_ROI,:,trial);
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
        [spindlePowerInBininTrial(trial,:),xedges_Placebo] = Phase_coup(FilteredEEG_SO, out_SO, out_SS, mode);
        
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
    v_MeanDensitySpindlesPlacebo(subj) = nanmean(v_DensitySpindles); %you can do either mean or sum
    v_MeanPowerSpindlesPlacebo(subj) = nanmean(v_MeanPowerSpindles); %you can do either mean or sum
    spindleInBininSubject_Placebo(subj,:) = nanmean(spindleInBininTrial,1); %you can do either mean or sum
    spindlePowerInBininSubject_Placebo(subj,:) = nanmean(spindlePowerInBininTrial,1); %you can do either mean or sum
end

%% Plots
figure
histogram(v_MeanDensitySpindlesOdor)
hold on
histogram(v_MeanDensitySpindlesPlacebo)
title('Density')


figure
histogram(v_MeanPowerSpindlesOdor)
hold on
histogram(v_MeanPowerSpindlesPlacebo)
title('Power')


%% Polar Plot for Density
meanSpindleInBin_Placebo = nanmean(spindleInBininSubject_Placebo,1);
meanSpindleInBin_Odor = nanmean(spindleInBininSubject_Odor,1);

spindleInBininSubject.Placebo = spindleInBininSubject_Placebo;
spindleInBininSubject.Odor = spindleInBininSubject_Odor;
save([cd, 'MeanSpindlesInBin'],'spindleInBininSubject');
figure;
polarhistogram('BinEdges', xedges_Placebo, 'BinCounts', meanSpindleInBin_Placebo,'FaceColor', 'magenta',...
    'FaceAlpha',.3);
thetaticks(0:90:315); %rlim([0 .02]);
pax = gca; pax.ThetaAxisUnits = 'radians';
pax.FontSize = 12; pax.GridColor = 'red';
title('SW Phase-Spindle Onset coupling - Placebo');

figure;
polarhistogram('BinEdges', xedges_Odor, 'BinCounts', meanSpindleInBin_Odor,'FaceColor', 'magenta',...
    'FaceAlpha',.3);
thetaticks(0:90:315); %rlim([0 .02]);
pax = gca; pax.ThetaAxisUnits = 'radians';
pax.FontSize = 12; pax.GridColor = 'red';
title('SW Phase-Spindle Onset coupling - Odor');

%% Polar Plot for Power

meanSpindlePowerInBin_Placebo = nanmean(spindlePowerInBininSubject_Placebo,1);
meanSpindlePowerInBin_Odor = nanmean(spindlePowerInBininSubject_Odor,1);

spindlePowerInBininSubject.Placebo = spindlePowerInBininSubject_Placebo;
spindlePowerInBininSubject.Odor = spindlePowerInBininSubject_Odor;
save([cd, 'MeanSpindlesPowerInBin'],'spindlePowerInBininSubject');
figure;
polarhistogram('BinEdges', xedges_Placebo, 'BinCounts', meanSpindlePowerInBin_Placebo,'FaceColor', 'magenta',...
    'FaceAlpha',.3);
thetaticks(0:90:315); rlim([0 2.2]);
pax = gca; pax.ThetaAxisUnits = 'radians';
pax.FontSize = 12; pax.GridColor = 'red';
title('SW Phase-Power Spindle Onset coupling - Placebo');

figure;
polarhistogram('BinEdges', xedges_Odor, 'BinCounts', meanSpindlePowerInBin_Odor,'FaceColor', 'magenta',...
    'FaceAlpha',.3);
thetaticks(0:90:315); rlim([0 2.2]);
pax = gca; pax.ThetaAxisUnits = 'radians';
pax.FontSize = 12; pax.GridColor = 'red';
title('SW Phase-Power Spindle Onset coupling - Odor');




%% Save useful information

if strcmp(ChanVsSource,'Channel')
    save(strcat('spindleInBin_',str_ChanSO,'vs',str_ChanSS,'_',OnVsOff), 'spindleInBinTotalOdor','spindleInBinTotalPlacebo');
    save(strcat('spindlePowerInBin',str_ChanSO,'vs',str_ChanSS,'_',OnVsOff), 'spindlePowerInBinTotalOdor','spindlePowerInBinTotalPlacebo');
    save(strcat('DensitySpindles',str_ChanSO,'vs',str_ChanSS,'_',OnVsOff), 'DensitySpindlesTotalOdor','DensitySpindlesTotalPlacebo');
    save(strcat('PowerSpindles',str_ChanSO,'vs',str_ChanSS,'_',OnVsOff), 'PowerSpindlesTotalOdor','PowerSpindlesTotalPlacebo');   
else  
    save(strcat('spindleInBin',str_ROI_SO,'vs',str_ROI_SS,'_',OnVsOff), 'spindleInBinTotalOdor','spindleInBinTotalPlacebo');
    save(strcat('spindlePowerInBin',str_ROI_SO,'vs',str_ROI_SS,'_',OnVsOff), 'spindlePowerInBinTotalOdor','spindlePowerInBinTotalPlacebo');
    save(strcat('DensitySpindles',str_ROI_SO,'vs',str_ROI_SS,'_',OnVsOff), 'DensitySpindlesTotalOdor','DensitySpindlesTotalPlacebo');
    save(strcat('PowerSpindles',str_ROI_SO,'vs',str_ROI_SS,'_',OnVsOff), 'PowerSpindlesTotalOdor','PowerSpindlesTotalPlacebo');  
end
