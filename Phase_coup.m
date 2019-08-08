function [spindleInBinTrial,xedges] = Phase_coup(FilteredEEG_ctx, out_ctx, out_th)

xhilb = hilbert(FilteredEEG_ctx);
xphase = angle(xhilb); 

spindles  = zeros(length(xphase),1); 
nedges = linspace(-pi,pi,15);

if isnan(out_th.trialinfo(:,2))
    spindleInBinTrial = nan(length(nedges)-1,1);
    xedges = nedges;
    return 
end

spindles(out_th.trialinfo(:,2)) = 1; 
spindles(spindles==1) = out_th.trialinfo(:,12);

for sw = 1:size(out_ctx.trialinfo,1)
phase_sw = xphase(out_ctx.trialinfo(sw,2):out_ctx.trialinfo(sw,4));
spindles_i = spindles(out_ctx.trialinfo(sw,2):out_ctx.trialinfo(sw,4));

[N,xedges,bin] = histcounts(phase_sw, nedges); 

for E=1:length(xedges)-1
    idxs = find(phase_sw>=xedges(E) & phase_sw<xedges(E+1));   
    spindleInBin(sw, E) = mean(spindles_i(idxs)); 
end
end 
spindleInBinTrial = sum(spindleInBin,1);

% polarhistogram('BinEdges', xedges, 'BinCounts', spindleInBinTrial,'FaceColor', 'magenta',...
%     'FaceAlpha',.3);
% thetaticks(0:90:315);
% pax = gca; pax.ThetaAxisUnits = 'radians';
% pax.FontSize = 12; pax.GridColor = 'red';

