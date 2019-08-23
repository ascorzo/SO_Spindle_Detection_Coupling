function [spindleInBinTrial,xedges] = Phase_coup(FilteredEEG_SO, out_SO, out_SS,mode)

% Calculation of Envelope and angle of the Filtered EEG
v_xhilb = hilbert(FilteredEEG_SO);
v_xphase = angle(v_xhilb); 

% Initialization of vector of spindles as nan, in order to do meannan of
% power in the future
v_spindles  = nan(length(v_xphase),1); 

% Division in bins of the 2pi phase of the SO, to detect the spindles
% coupled
% v_nedges = linspace(-pi,pi,15); % Change here Bins ***
v_nedges = linspace(-pi,pi,9);

%If there are no spindles, then steps out the function and returns a vector
%of nan
if isnan(out_SS.trialinfo(:,2))
    spindleInBinTrial = nan(length(v_nedges)-1,1);
    xedges = v_nedges;
    return 
end

% Make the v_spindles vector equal to ones when there is a spindle onset
v_spindles(out_SS.trialinfo(:,2)) = 1; 

% If we want to check the power, then, multiplies the v_spindle vector by
% the power of each spindle associated onset
if strcmp(mode,'Power') == 1
    v_spindles(v_spindles==1) = out_SS.trialinfo(:,12);
end


for sw = 1:size(out_SO.trialinfo,1)
    % calculates the phase in which each SO falls  
    %v_xphase(startpointSO:endpointSO)
    phase_sw = v_xphase(out_SO.trialinfo(sw,2):out_SO.trialinfo(sw,4));
    
    %Calculates the spindles onset or power located in each SO
    spindles_i = v_spindles(out_SO.trialinfo(sw,2):out_SO.trialinfo(sw,4));
    
    
    [N,xedges,bin] = histcounts(phase_sw, v_nedges);
    
    for Edge=1:length(xedges)-1
        %Find indexes of each phase bin
        idxs = find(phase_sw>=xedges(Edge) & phase_sw<xedges(Edge+1));
        if strcmp(mode,'Power') == 1
            %Calculates the mean of spindles in the indexes of each phase bin
            spindleInBin(sw, Edge) = nanmean(spindles_i(idxs));
        else
            %Calculates the total of spindles in the indexes of each phase bin
            spindleInBin(sw, Edge) = nansum(spindles_i(idxs));
        end
    end
end
if strcmp(mode,'Power') == 1
    % calculates the mean of spindles in each bin for all the SO detected
    spindleInBinTrial = nanmean(spindleInBin,1);
else
    % calculates the total of spindles in each bin for all the SO detected
    spindleInBinTrial = nansum(spindleInBin,1);
end

