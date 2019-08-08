%% Define some important variables
s_TotalTimeSec = 30; %total time in seconds we are analysing
s_fs = 1000; % sample frequency
s_TimeSam = s_TotalTimeSec * s_fs; %total time in samples

% In case the time is unbalanced before and after the stimuli, specify the times.
s_TimeAfterCero = 15;
s_TimeBeforeCero = 15;

str_Chan = "C3";
%Region of interest for Cortex:
%Possible options: {'Frontal_Inf_Oper L';'Frontal_Inf_Oper R';'Frontal_Inf_Orb_2 L';'Frontal_Inf_Orb_2 R';'Frontal_Inf_Tri L';'Frontal_Inf_Tri R';'Frontal_Med_Orb L';'Frontal_Med_Orb R';'Frontal_Mid_2 L';'Frontal_Mid_2 R';'Frontal_Sup_2 L';'Frontal_Sup_2 R';'Frontal_Sup_Medial L';'Frontal_Sup_Medial R'}
str_ROI_ctx = 'Frontal_Sup_Medial L';

%Region of interest for Thalamus:
%Possible options: {'Thalamus L';'Thalamus R'}
str_ROI_th = 'Thalamus L';


%Other possible options of regions to compare:h
%{'Amygdala L';'Amygdala R';'Angular L';'Angular R';'Calcarine L';'Calcarine R';
%'Caudate L';'Caudate R';'Cerebelum_3 R';'Cerebelum_4_5 L';'Cerebelum_4_5 R';
%'Cerebelum_6 L';'Cerebelum_6 R';'Cerebelum_Crus1 L';'Cerebelum_Crus1 R';
%'Cerebelum_Crus2 L';'Cingulate_Ant L';'Cingulate_Ant R';'Cingulate_Mid L';
%'Cingulate_Mid R';'Cingulate_Post L';'Cingulate_Post R';'Cuneus L';'Cuneus R';
%'Fusiform L';'Fusiform R';'Heschl L';'Heschl R';'Hippocampus L';'Hippocampus R';
%'Insula L';'Insula R';'Lingual L';'Lingual R';'Occipital_Inf L';'Occipital_Inf R';
%'Occipital_Mid L';'Occipital_Mid R';'Occipital_Sup L';'Occipital_Sup R';'OFCant L';
%'OFCant R';'OFClat L';'OFClat R';'OFCmed L';'OFCmed R';'OFCpost L';'OFCpost R';
%'Olfactory L';'Olfactory R';'Pallidum L';'Pallidum R';'Paracentral_Lobule L';
%'Paracentral_Lobule R';'ParaHippocampal L';'ParaHippocampal R';'Parietal_Inf L';
%'Parietal_Inf R';'Parietal_Sup L';'Parietal_Sup R';'Postcentral L';'Postcentral R';
%'Precentral L';'Precentral R';'Precuneus L';'Precuneus R';'Putamen L';'Putamen R';
%'Rectus L';'Rectus R';'Rolandic_Oper L';'Rolandic_Oper R';'Supp_Motor_Area L';
%'Supp_Motor_Area R';'SupraMarginal L';'SupraMarginal R';'Temporal_Inf L';
%'Temporal_Inf R';'Temporal_Mid L';'Temporal_Mid R';'Temporal_Pole_Mid L';
%'Temporal_Pole_Mid R';'Temporal_Pole_Sup L';'Temporal_Pole_Sup R';'Temporal_Sup L';'Temporal_Sup R'}


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
                findChannel = find(strcmp(fileChan.Channel.label, str_Chan));
                DataOdorTrial_ctx = DataOdorChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_ctx);
                DataOdorTrial_ctx = DataOdorSource(s_ROI,:,trial);
        end
        
        Label = 'ctx';
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
        ylim([-.3 .3])
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Spindle detection in thalamus
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Store choosen region for later saving
        LabelCortex = Label;
        
        switch ChanVsSource
            case 'Channel'
                findChannel = find(strcmp(fileChan.Channel.label, str_Chan));
                DataOdorTrial_th = DataOdorChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_th);
                DataOdorTrial_th = DataOdorSource(s_ROI,:,trial);
        end
        Label = 'Thalamus';
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
        ylim([-.3 .3])
        %---------------------------------------------------
        
        [spindleInBininTrial(trial,:),xedges_Odor] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th);
        v_DensitySpindles(trial) = size(out_th.trialinfo,1);%/s_TimeAfterCero; % Play with density or just number (Density as # of spindles in time)
    end
    v_MeanDensitySpindlesOdor(subj) = nanmean(v_DensitySpindles); %you can do either mean or sum
    spindleInBininSubject_Odor(subj,:) = nanmean(spindleInBininTrial,1); %you can do either mean or sum
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
                findChannel = find(strcmp(fileChan.Channel.label, str_Chan));
                DataPlaceboTrial_ctx = DataPlaceboChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_ctx);
                DataPlaceboTrial_ctx = DataPlaceboSource(s_ROI,:,trial);
        end
        
        Label = 'ctx';
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
                findChannel = find(strcmp(fileChan.Channel.label, str_Chan));
                DataPlaceboTrial_th = DataPlaceboChan(findChannel,:,trial);
            case 'Source'
                s_ROI = strcmp({fileSource.Atlas.Scouts.Label},str_ROI_th);
                DataPlaceboTrial_th = DataPlaceboSource(s_ROI,:,trial);
        end
    
        Label = 'Thalamus';
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
        
        [spindleInBininTrial(trial,:),xedges_Placebo] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th);
        v_DensitySpindles(trial) = size(out_th.trialinfo,1);%/s_TimeAfterCero; % Play with density or just number (Density as # of spindles in time)
    end
    v_MeanDensitySpindlesPlacebo(subj) = nansum(v_DensitySpindles); %you can do either mean or sum
    spindleInBininSubject_Placebo(subj,:) = nanmean(spindleInBininTrial,1); %you can do either mean or sum
end

%% Plots
figure
histogram(v_MeanDensitySpindlesOdor,10)
hold on
histogram(v_MeanDensitySpindlesPlacebo,10)

figure
histogram((v_MeanDensitySpindlesOdor-v_MeanDensitySpindlesPlacebo),22);

%%
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