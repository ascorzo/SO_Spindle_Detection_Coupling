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

ChanVsSource = questdlg('Import channel recording data or source-computed cortex data?', ...
    'Channel or Source?', ...
    'Channel','Source','Channel');

OnVsOff = questdlg('Select if you want to evaluate off condition or odor condition', ...
    'Odor On or Odor Off?', ...
    'On','Off','On');


%%

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
            
            %if numel(FilesListChanOdor) ~= numel(FilesListSourceOdor) || numel(FilesListChanPlacebo) ~= numel(FilesListSourcePlacebo)
            %    error('Mismatch in numbers of datasets')
            %end
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
    %otherwise
     %   FilesListSourceOdor = dir([pathNameSource,'*Odor.mat']);
      %  FilesListSourcePlacebo = dir([pathNameSource,'*Placebo.mat']);
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
    [indx,tf] = listdlg('PromptString','Select a source: for Slow Oscillations',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ChanSO = Sources(indx);
    
    [indx,tf] = listdlg('PromptString','Select a source: for Sleep Spindles',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ChanSS = Sources(indx);
    
else
    FileTemp = subjectFiles{1};
    Sources = string({FileTemp.Atlas.Scouts.Label});
    [indx,tf] = listdlg('PromptString','Select a source for SO:',...
        'SelectionMode','single',...
        'ListString',Sources);
    str_ROI_SO = Sources(indx);
    [indx,tf] = listdlg('PromptString','Select a source for Spindles:',...
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
        %  SW detection in CTX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ChanVsSource
            case 'Channel'
                findChannel = find(strcmp(fileChan.Channel.label, str_ChanSO));
                DataOdorTrial_ctx = DataOdorChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                DataOdorTrial_ctx = DataOdorSource(s_ROI,:,trial);
        end
        
        Label = 'SO';
        SleepScoring = ones(length(DataOdorTrial_ctx),1);
        ScoringArtefacts = zeros(length(DataOdorTrial_ctx),1);
        if strcmp(OnVsOff,'On') == 1
            ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
        else
            ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
        end
        
        [FilteredEEG_ctx, out_ctx] = SODetection(DataOdorTrial_ctx', SleepScoring, ScoringArtefacts, Label, s_fs);
        
        v_time = (0:1/s_fs:(length(FilteredEEG_ctx)-1)/s_fs) - s_TimeBeforeCero;

        %------------------------------------------------------------------
        % Plot in case you want to check the detection
        %------------------------------------------------------------------
        subplot(2,1,1);
        plot(v_time,FilteredEEG_ctx);
        for i=1:size(out_ctx.trialinfo,1)
            line([out_ctx.trialinfo(i,2)/out_ctx.fsample out_ctx.trialinfo(i,4)/out_ctx.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
        end
        %ylim([-.3 .3])
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Spindle detection in thalamus
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Store choosen region for later saving
        LabelCortex = Label;
        
        switch ChanVsSource
            case 'Channel'
                findChannel = find(strcmp(fileChan.Channel.label, str_ChanSS));
                DataOdorTrial_th = DataOdorChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                DataOdorTrial_th = DataOdorSource(s_ROI,:,trial);
        end
        Label = 'SS';
        SleepScoring = ones(length(DataOdorTrial_th),1);
        ScoringArtefacts = zeros(length(DataOdorTrial_th),1);
        
        if strcmp(OnVsOff,'On') == 1
            ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
        else
            ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
        end
        
        [FilteredEEG_th, out_th] = SpindleDetection(DataOdorTrial_th', SleepScoring, ScoringArtefacts, Label, s_fs);
        
        %---------------------------------------------
        subplot(2,1,2);
        plot(v_time,FilteredEEG_th);
        for i=1:size(out_th.trialinfo,1)
            line([out_th.trialinfo(i,2)/out_th.fsample out_th.trialinfo(i,4)/out_th.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
        end
        %ylim([-.3 .3])
        %---------------------------------------------------
                
        mode = 'Density';
        [spindleInBininTrial(trial,:),xedges_dor] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th, mode);
        
        mode = 'Power';
        [spindlePowerInBininTrial(trial,:),xedges_Odor] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th, mode);
        
        v_DensitySpindles(trial) = size(out_th.trialinfo,1);
        v_MeanPowerSpindles(trial) = mean(out_th.trialinfo(:,12));
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
        %  SW detection in CTX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ChanVsSource
            case 'Channel'
                findChannel = find(strcmp(fileChan.Channel.label, str_ChanSO));
                DataPlaceboTrial_ctx = DataPlaceboChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SO);
                DataPlaceboTrial_ctx = DataPlaceboSource(s_ROI,:,trial);
        end
        
        Label = 'SO';
        SleepScoring = ones(length(DataPlaceboTrial_ctx),1);
        ScoringArtefacts = zeros(length(DataPlaceboTrial_ctx),1);
        if strcmp(OnVsOff,'On') == 1
            ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
        else
            ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
        end
        
        [FilteredEEG_ctx, out_ctx] = SODetection(DataPlaceboTrial_ctx', SleepScoring, ScoringArtefacts, Label, s_fs);
        
        v_time = (0:1/s_fs:(length(FilteredEEG_ctx)-1)/s_fs) - s_TimeBeforeCero;
        
        %------------------------------------------------------------------
        % Plot in case you want to check the detection
        %------------------------------------------------------------------
        
                subplot(2,1,1);
                plot(v_time,FilteredEEG_ctx);
                for i=1:size(out_ctx.trialinfo,1)
                    line([out_ctx.trialinfo(i,2)/out_ctx.fsample out_ctx.trialinfo(i,4)/out_ctx.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
                end
                ylim([-.3 .3])
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Spindle detection in thalamus
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Store choosen region for later saving
        LabelCortex = Label;
        
        switch ChanVsSource
            case 'Channel'
                findChannel = find(strcmp(fileChan.Channel.label, str_ChanSS));
                DataPlaceboTrial_th = DataPlaceboChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_SS);
                DataPlaceboTrial_th = DataPlaceboSource(s_ROI,:,trial);
        end
    
        Label = 'SS';
        SleepScoring = ones(length(DataPlaceboTrial_th),1);
        ScoringArtefacts = zeros(length(DataPlaceboTrial_th),1);
        if strcmp(OnVsOff,'On') == 1
            ScoringArtefacts(1:s_TimeBeforeCero*s_fs) = 1;
        else
            ScoringArtefacts(s_TimeBeforeCero*s_fs:end) = 1;
        end
        
        [FilteredEEG_th, out_th] = SpindleDetection(DataPlaceboTrial_th', SleepScoring, ScoringArtefacts, Label, s_fs);
        
        %----------------------------------------------
        subplot(2,1,2);
        plot(v_time,FilteredEEG_th);
        for i=1:size(out_th.trialinfo,1)
            line([out_th.trialinfo(i,2)/out_th.fsample out_th.trialinfo(i,4)/out_th.fsample]-s_TimeBeforeCero,[-.2,-.2], 'Color','red')
        end
        ylim([-.3 .3])
        %-----------------------------------------
        mode = 'Density';
        [spindleInBininTrial(trial,:),xedges_Placebo] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th, mode);
        
        mode = 'Power';
        [spindlePowerInBininTrial(trial,:),xedges_Placebo] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th, mode);
        
        v_DensitySpindles(trial) = size(out_th.trialinfo,1);%/s_TimeAfterCero; % Play with density or just number (Density as # of spindles in time)
        v_MeanPowerSpindles(trial) = mean(out_th.trialinfo(:,12));
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
Spindles.Run = datetime;
Spindles.ColumnNames = {'Odor', 'Placebo'};
Spindles.subjPwrOdor = spindleInBininSubject_Odor;
Spindles.subjPwrPlacebo = spindleInBininSubject_Placebo;
Spindles.MeanPwr = [meanSpindleInBin_Odor' meanSpindleInBin_Placebo'];
% Need here the number computation from the old script, please
%Spindles.subjNumberOdor = ;
%Spindles.subjNumberPlacebo = ;
%Spindles.MeanNumber = ;
Spindles.DatasetInfo.s_TotalTimeSec = s_TotalTimeSec;
Spindles.DatasetInfo.Samplerate = s_fs;
Spindles.DatasetInfo.NumberSubjects = [numel(FilesListSourceOdor), numel(FilesListSourcePlacebo)];
Spindles.DatasetInfo.EpochSize = [s_TimeBeforeCero, s_TimeAfterCero];
Spindles.ROI = {LabelCortex, LabelSubCtx};

save([cd slashSys 'run_SODetection'], 'Spindles');