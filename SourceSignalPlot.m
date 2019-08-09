clear all

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

%%

if ~exist('pathNameChan','var')
    % path to datasets containing channel data
    pathNameChan = [uigetdir(cd,'Choose the folder that contains the CHANNEL datasets'), slashSys];
    addpath(pathNameChan)
    FilesListChanOdor = dir([pathNameChan,'*Odor.mat']);
    FilesListChanPlacebo = dir([pathNameChan,'*Placebo.mat']);  
    
end

if ~exist('pathNameSource','var')
    
    %path to datasets with source estimation
    pathNameSource = [uigetdir(cd,'Choose the folder that contains the SOURCE datasets'), slashSys];
    addpath(pathNameSource)
    FilesListSourceOdor = dir([pathNameSource,'*Odor.mat']);
    FilesListSourcePlacebo = dir([pathNameSource,'*Placebo.mat']);
end
            
if numel(FilesListChanOdor) ~= numel(FilesListSourceOdor) || numel(FilesListChanPlacebo) ~= numel(FilesListSourcePlacebo)
   error('Mismatch in numbers of datasets')
end

%% Select the signals you want to plot

FileTemp = load(FilesListChanOdor(1).name);
Sources = FileTemp.Channel.label;
[indx,tf] = listdlg('PromptString','Select a channnel for scalp signal',...
    'SelectionMode','single','Listsize',[300,300],...
    'ListString',Sources);
str_scalp = Sources(indx);

FileTemp = load([pathNameSource FilesListSourceOdor(1).name]);
Sources = string({FileTemp.Atlas.Scouts.Label});
[indx,tf] = listdlg('PromptString','Select a source for cortex signal',...
    'SelectionMode','single','Listsize',[300,300],...
    'ListString',Sources);
str_ctx = Sources(indx);

[indx,tf] = listdlg('PromptString','Select a source for subcortical signal',...
    'SelectionMode','single','Listsize',[300,300],...
    'ListString',Sources);
str_subct = Sources(indx);


%% Loading Odor sets into memory
for Subj = 1:numel(FilesListChanOdor)
    subjectFilesChan{Subj,1} = load([pathNameChan FilesListChanOdor(Subj).name]);
    subjectFilesSource{Subj,1} = load([pathNameSource FilesListSourceOdor(Subj).name]);
    
    fileChan = load(FilesListChanOdor(Subj).name);
    s_NumChannels = length(fileChan.Channel.number);
    v_Time = fileChan.Channel.times;
    s_Trials = fileChan.Channel.trials;
    s_NumScouts = length(subjectFilesSource{Subj,1}.Atlas.Scouts);
    
    DataOdorChan = double(fileChan.Channel.data);
    DataOdorSource = reshape(subjectFilesSource{Subj,1}.Value,[s_NumScouts,s_TimeSam,s_Trials]);
    
    findChannel = find(strcmp(fileChan.Channel.label, str_scalp));
    s_ROI_ctx = strcmp({subjectFilesSource{Subj,1}.Atlas.Scouts.Label},str_ctx);
    s_ROI_subct = strcmp({subjectFilesSource{Subj,1}.Atlas.Scouts.Label},str_subct);
    
    
    for trial = 1:s_Trials
        DataOdorTrial_scalp = DataOdorChan(findChannel,:,trial);
        DataOdorTrial_ctx = DataOdorSource(s_ROI_ctx,:,trial);
        DataOdorTrial_subct = DataOdorSource(s_ROI_subct,:,trial);
        
        figure
        subplot(3,1,1)
        plot(v_Time,DataOdorTrial_scalp)
        
        ymin = min([DataOdorTrial_ctx,DataOdorTrial_subct]);
        ymax = max([DataOdorTrial_ctx,DataOdorTrial_subct]);
        
        subplot(3,1,2)
        plot(v_Time,DataOdorTrial_ctx)
        ylim([ymin,ymax])
        
        subplot(3,1,3)
        plot(v_Time,DataOdorTrial_subct)
        ylim([ymin,ymax])        
        
        %% Spindle detection
        Label = 'Spindles';
        SleepScoring = ones(length(DataOdorTrial_ctx),1);
        ScoringArtefacts = zeros(length(DataOdorTrial_ctx),1);
        
        
        [FilteredEEG_scalp, out_scalp] = SpindleDetection(DataOdorTrial_scalp', SleepScoring, ScoringArtefacts, Label, s_fs);
        v_time = (0:1/s_fs:(length(FilteredEEG_scalp)-1)/s_fs) - s_TimeBeforeCero;
        
        [FilteredEEG_ctx, out_ctx] = SpindleDetection(DataOdorTrial_ctx', SleepScoring, ScoringArtefacts, Label, s_fs);
        [FilteredEEG_subct, out_subct] = SpindleDetection(DataOdorTrial_subct', SleepScoring, ScoringArtefacts, Label, s_fs);

        %------------------------------------------------------------------
        % Plot in case you want to check the detection
        %------------------------------------------------------------------
        figure
        subplot(3,1,1)
        plot(v_time,FilteredEEG_scalp)
        
        ymin = min([FilteredEEG_ctx;FilteredEEG_subct]);
        ymax = max([FilteredEEG_ctx;FilteredEEG_subct]);
        
        subplot(3,1,2)
        plot(v_time,FilteredEEG_ctx)
        ylim([ymin,ymax])
        
        subplot(3,1,3);
        plot(v_time,FilteredEEG_subct)
        ylim([ymin,ymax])
        
    end 
end
