function [FilteredEEG, out] = SpindleDetection(RecordingValues, SleepScoring, ScoringArtefacts, Label, fsample)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spindle Detection Algorith for EEG/LFP siganls
% INPUT:    RecordingValues: Matrix containing Raw EEG or LFP values
%           SleepScoring: Matrix containing Sleep Scoring for each data point
%           ScoringArtefact: Matrix coding the artefacts (1) of the signal
%           for each data point.
%           Label: CHAR with the name of the Channel used.
%           fsample: sampling rate.
% OUT:      Output contains detection info in the following way:
%           out.label           = CELL containg the name of the Channel(s)
%           out.fsample         = sampling rate
%           out.trial           = Event(s) raw signal
%           out.time            = time/length of the Trial(s)
%           out.trialinfo       = Detailed info about each one of the trials.
%               out.trialinfo(:,1) = To which time bin corresponds that particular event 
%               out.trialinfo(:,2) = startTime: begin of event 
%               out.trialinfo(:,3) = midTime: center of event
%               out.trialinfo(:,4) = endTime: spindle: negative zero crossing
%               out.trialinfo(:,5) = duration: duration from start to end in seconds (spindle: betweet the two threshild crossings)
%               out.trialinfo(:,6) = maxTime: time of maximum (spindle: largest negative trough during spindle) in datapoints of original dataset
%               out.trialinfo(:,7) = minTime: time of minimum (spindle: largest positive peak during spindle) in datapoints of original dataset
%               out.trialinfo(:,8) = minAmp: amplitude of maximum (spindle: largest negative trough during spindle) in µV
%               out.trialinfo(:,9) = maxAmp: amplitude of minimum (spindle: largest positive peak during spindle) in µV
%               out.trialinfo(:,10)= p2pAmp: peak-to-peak amplitude (spindle: largest peak to largest trough) 
%               out.trialinfo(:,11)= p2pTime: time in seconds from peak-to-peak (spindle: abs(min-max)) 
%               out.trialinfo(:,12)= Power
%               out.trialinfo(:,13)= Frequency
%           out.vector          = Vector containing the presence of the events
%           out.sleeptimeinsec  = Time in SECONDS used for the Detections (excl. time spent in Artefacts)
%
% Authors:  Carlos N. Oyanedel - Uni Tübingen, Germany
%           Niels Niethard - Uni Tübingen, Germany
%           Thanks to Dr. Hong-Viet Ngo, University of Birmingham, UK
% Contact:  jan.born@uni-tuebingen.de
%
% Date:     17.09.2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic Parameters
InvalidSpindles = 0;

% Spindle Duration in Seconds
MinSpindleDur = 0.4;
MaxSpindleDur = 3;

% Filters in Hz
LowPassFilter       = 16;
HighPassFilter      = 7;

% Sleep Stage of Interest
%1: WAKE 2: NREM; 3:REM; 4:PREM 
SleepStages     = 1:4;
SleepStageToUse = 1;

% STD above mean
StdValue        = .2; %The amount of STD above the mean for the detection

% Trial Window in seconds
TrialWindow  = 2.5; % Window before and after the Min Peak

%% Filtering the signal and calculating the envelope
[d,e]               = butter(3,2*LowPassFilter/fsample,'low'); %Use butterworth filter 3rd order. 15Hz cutoff
FilteredEEG         = filtfilt(d,e,RecordingValues); %Filter Signal
[d,e]               = butter(3,2*HighPassFilter/fsample,'high'); %Use butterworth filter 3rd order. 8Hz cutoff
FilteredEEG         = filtfilt(d,e,FilteredEEG);
FilteredEEGEnvelope = abs(hilbert(FilteredEEG));

%% Smoothing the FilteredEEGEnvelope
run SmoothFilteredEEGEnvelope.m
%% Mean and STD of the SmoothedFilteredEEGEnvelope of an specific Sleep Stage
ValidSleepScoring                       = SleepScoring;
ValidSleepScoring(ScoringArtefacts==1)  = 0; %Delete data with artefacts
SmoothedFilteredEEGEnvelopeMean         = mean(SmoothedFilteredEEGEnvelope(find(ValidSleepScoring==SleepStageToUse),1));
SmoothedFilteredEEGEnvelopeSTD          = std(SmoothedFilteredEEGEnvelope(find(ValidSleepScoring==SleepStageToUse)));

%% SpindleVector.
SpindleVector   = zeros(length(RecordingValues),1);

for i=1:length(RecordingValues); %Detecting Spindles crossing the power threshold
    if SmoothedFilteredEEGEnvelope(i) > (SmoothedFilteredEEGEnvelopeMean+StdValue*SmoothedFilteredEEGEnvelopeSTD)
        SpindleVector(i) = 1;
    end
end

%% If want to delete event detected in different sleep stages
% Delete detected values in the other Sleep Stages
SpindleVector(1,1)                      = 0;
SleepStagetoRemove                      = SleepStages(SleepStages~=SleepStageToUse);
SpindleVector(ValidSleepScoring==SleepStagetoRemove(1)) = 0; %
SpindleVector(ValidSleepScoring==SleepStagetoRemove(2)) = 0; %
SpindleVector(ValidSleepScoring==SleepStagetoRemove(3)) = 0; %
SpindleVector(ValidSleepScoring==8)     = 0; %In Artefacts 8
SpindleVector(ScoringArtefacts==1)      = 0; %In Artefacts
SpindleVector(ValidSleepScoring==0)     = 0; %In case of mistmatch in extending the Sleep Scoring ,i.e. DiffScoringRec>0

%% Detecting Event Beg and End
% Position of values detected
SpindleVectorLoc    = find(SpindleVector==1); 

% Beg and End of each event.
SpindleEnd =[];
SpindleBeg =[];

for i=2:length(SpindleVectorLoc)-1
    if SpindleVectorLoc(i) - SpindleVectorLoc(i-1) > 1
        SpindleBeg = [SpindleBeg,SpindleVectorLoc(i)];
    end
    if SpindleVectorLoc(i+1) - SpindleVectorLoc(i) > 1 
        SpindleEnd = [SpindleEnd,SpindleVectorLoc(i)];
    end
end

if numel(SpindleBeg) == numel(SpindleEnd)
    SpindleBeg          = [SpindleVectorLoc(1), SpindleBeg];
    SpindleEnd          = [SpindleEnd, SpindleVectorLoc(end)];
    
elseif numel(SpindleBeg) < numel(SpindleEnd)
    SpindleBeg          = [SpindleVectorLoc(1), SpindleBeg];
    
elseif numel(SpindleBeg) > numel(SpindleEnd)
    SpindleEnd          = [SpindleEnd, SpindleVectorLoc(end)];
end
SpindlesDetected    = [SpindleBeg; SpindleEnd];
%% Duration Threshold
% Now look for spindles with a Min and Max length
SpindleDuration = diff(SpindlesDetected);

MinSpindleDur = MinSpindleDur*fsample;
MaxSpindleDur = MaxSpindleDur*fsample;

j = 1;
for i = 1:size(SpindlesDetected,2)
    if SpindleDuration(i) <= MaxSpindleDur && SpindleDuration(i) >= MinSpindleDur %Checking which ones fulfill the duration criteria
        ValidSpindles(1,j)  = SpindlesDetected(1,i);
        ValidSpindles(2,j)  = SpindlesDetected(2,i);
        j                   = j+1;
    end
end
if ~exist("ValidSpindles")
    ValidSpindles = SpindlesDetected;
    InvalidSpindles = 1;
end
    
ValidSpindleDurationSec = (diff(ValidSpindles))./fsample;
% Peaks
PeakSpindleNumber = zeros(1,size(ValidSpindles,2));
for i = 1:size(ValidSpindles,2);
    PeakSpindleNumber(:,i) = length(findpeaks(FilteredEEG(ValidSpindles(1,i):ValidSpindles(2,i))));
end
% Frequency
tmpFreq = PeakSpindleNumber./ValidSpindleDurationSec;


%% Valid Ripple Vector
ValidSpindleVector = zeros(1,size(FilteredEEG,1));
for i=1:size(ValidSpindles,2)
    ValidSpindleVector(1,ValidSpindles(1,i):ValidSpindles(2,i))=1;
end

%% Spindle Power
tmppower = zeros(size(ValidSpindles,2),1);
for iPow = 1:size(ValidSpindles,2);
    tmppower(iPow,1) = trapz(SmoothedFilteredEEGEnvelope(ValidSpindles(1,iPow):ValidSpindles(2,iPow)));
end

%% Detecting Max Values and position for the Envelope of each Valid Spindle
MaxSpindleEnvelopeVal = []; % Local Maxima Value Matrx for the Envelope
MaxSpindleEnvelopeLoc = []; % Local Maxima Position for the Envelope
for iMax = 1:size(ValidSpindles,2);
    [MaxSpindleEnvelopeVal(:,iMax), MaxSpindleEnvelopeLoc(:,iMax)] = max(SmoothedFilteredEEGEnvelope(ValidSpindles(1,iMax):ValidSpindles(2,iMax)));
end

%Positions related to the whole recording for the Envelope Valid Spindle Max Detected
for iPosRel = 1:size(ValidSpindles,2);
    ValidSpindleEnvelopeMaxLoc(iPosRel) = (ValidSpindles(1,iPosRel)+MaxSpindleEnvelopeLoc(iPosRel));
end

%% Detecting Local Max and Min
% Local Minima
MinValidSpindleVal = []; % Local Minima Value Matrx
MinValidSpindleLoc = []; % Local Minima Position
for iLocMin = 1:size(ValidSpindles,2);
    [MinValidSpindleVal(:,iLocMin), MinValidSpindleLoc(:,iLocMin)] = min(FilteredEEG(ValidSpindles(1,iLocMin):ValidSpindles(2,iLocMin)));
end

% Local Maxima
MaxValidSpindleVal = []; % Local Maxima Value Matrx
MaxValidSpindleLoc = []; % Local Maxima Position
for iLocMax = 1:size(ValidSpindles,2);
    [MaxValidSpindleVal(:,iLocMax), MaxValidSpindleLoc(:,iLocMax)] = max(FilteredEEG(ValidSpindles(1,iLocMax):ValidSpindles(2,iLocMax)));
end

% Positions related to the whole recording for the Min and Max for each Valid Spindle Detected
for iPosMinMax = 1:size(ValidSpindles,2);
    ValidSpindleMinLoc(iPosMinMax)  = (ValidSpindles(1,iPosMinMax) + MinValidSpindleLoc(iPosMinMax));
    ValidSpindleMaxLoc(iPosMinMax)  = (ValidSpindles(1,iPosMinMax) + MaxValidSpindleLoc(iPosMinMax));
end

% Valid Max2Min info (Time Difference)
Max2MinValidLoc     = [MinValidSpindleLoc; MaxValidSpindleLoc];
Max2MinValidDiffLoc = diff(Max2MinValidLoc);

% Valid Max2Min info (Amplitud)
Max2MinValidSpindles        = [MinValidSpindleVal; MaxValidSpindleVal];
Max2MinValidSpindlesDiff    = diff(Max2MinValidSpindles);

%% Grand Average
SpindleGrandAverage = zeros(size(ValidSpindles,2),5000);
out.trial           = {};
out.time            = {};
WindowGrandAverage  = TrialWindow*fsample; % Window before and after the Min Peak
Cut                 = 0; %In case the Trial Window for the last detected event is longer than the recording 
CutBeg              = 0; %In case the Trial Window for the first detected starts before the recording file

for i = 1:size(ValidSpindles,2);
    if ((ValidSpindles(1,i)+MinValidSpindleLoc(i))+WindowGrandAverage-1) > size(RecordingValues,1)
        Cut = 1;
    elseif ((ValidSpindles(1,i)+MinValidSpindleLoc(i))-WindowGrandAverage-1) < 0
        CutBeg = 1;        
    else
        SpindleGrandAverage(i,:)    = RecordingValues((ValidSpindles(1,i)+MinValidSpindleLoc(i))-WindowGrandAverage:(ValidSpindles(1,i)+MinValidSpindleLoc(i))+WindowGrandAverage-1);
        out.trial{1,i}              = SpindleGrandAverage(i,:)*1000;
        out.time{1,i}               = -TrialWindow:1/fsample:TrialWindow-1/fsample;
    end
end

%% Creating eventdata - Output
out.label       = {Label};
out.fsample     = fsample;
tmptrialinfo    = zeros(size(ValidSpindles,2),13);
for iTrial=1:size(ValidSpindles,2);
    tmptrialinfo(iTrial,1)   = (ceil(ValidSpindles(1,iTrial)/1800000)*30); %To which timebin corresponds
    tmptrialinfo(iTrial,2)   = ValidSpindles(1,iTrial); %startTime: begin of event (spindle: positive threshold crossing)
    tmptrialinfo(iTrial,3)   = (ValidSpindles(1,iTrial)+ValidSpindles(2,iTrial))/2; %midTime: center of event (spindle: largest trough)
    tmptrialinfo(iTrial,4)   = ValidSpindles(2,iTrial); %endTime: (spindle: negative zero crossing) 
    tmptrialinfo(iTrial,5)   = ValidSpindleDurationSec(iTrial); %duration: duration from start to end in seconds (spindle: betweet the two threshild crossings)
    tmptrialinfo(iTrial,6)   = ValidSpindleMinLoc(iTrial); %maxTime: time of maximum (spindle: largest negative trough during spindle) in datapoints of original dataset
    tmptrialinfo(iTrial,7)   = ValidSpindleMaxLoc(iTrial); %minTime: time of minimum (spindle: largest positive peak during spindle) in datapoints of original dataset
    tmptrialinfo(iTrial,8)   = MinValidSpindleVal(iTrial)*1000; %minAmp: amplitude of maximum (spindle: largest negative trough during spindle) in ï¿½V
    tmptrialinfo(iTrial,9)   = MaxValidSpindleVal(iTrial)*1000;%maxAmp: amplitude of minimum (spindle: largest positive peak during spindle) in ï¿½V
    tmptrialinfo(iTrial,10)  = Max2MinValidSpindlesDiff(iTrial)*1000; %p2pAmp: peak-to-peak amplitude (spindle: largest peak to largest trough)
    tmptrialinfo(iTrial,11)  = abs(Max2MinValidDiffLoc(iTrial))/fsample; %p2pTime: time in seconds from peak-to-peak (spindle: abs(min-max)) 
end
tmptrialinfo(:,12)  = tmppower'; %Power of each event detected
tmptrialinfo(:,13)  = tmpFreq'; %Freq of each event detected
if InvalidSpindles == 1
    tmptrialinfo = NaN(size(tmptrialinfo));
end
if Cut == 1
   tmptrialinfo = tmptrialinfo(1:size(out.time,2),:);
end

if CutBeg == 1
    out.trial   = out.trial(2:end);
    out.time    = out.time(2:end);
end

out.trialinfo       = tmptrialinfo;
out.vector          = ValidSpindleVector;
out.sleeptimeinsec  = ((size(find(SleepScoring == SleepStageToUse),1)) - (size(find(SleepScoring == SleepStageToUse & ScoringArtefacts == 1),1)))/fsample;

end