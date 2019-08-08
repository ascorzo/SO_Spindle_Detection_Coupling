% Smoothing FilteredEEGEnvelope signal.
% MovingAverages = input('Point window for smoothing the signal ');%Point window for smoothing the signal
MovingAverages = 200;
SmoothedFilteredEEGEnvelope = zeros(size(FilteredEEGEnvelope,1),1);

for i = 1:size(FilteredEEGEnvelope,1)
    if i < MovingAverages/2 % SmoothedMean for the first 100 data points
        X = FilteredEEGEnvelope(1:(round(MovingAverages/2)+(i-1)));

        MeanData = mean (X);
        SmoothedFilteredEEGEnvelope (i) = MeanData;
    end
    
    if i >= MovingAverages/2 && i < size(FilteredEEGEnvelope,1) - MovingAverages/2
        X = FilteredEEGEnvelope((i+1)-round(MovingAverages/2):round(MovingAverages/2)+(i-1));
 
        MeanData = mean (X);
        SmoothedFilteredEEGEnvelope (i) = MeanData;
    end
    if i >= size(FilteredEEGEnvelope,1) - MovingAverages/2 % SmoothedMean for the last 100 data points
        X = FilteredEEGEnvelope((i+1)-round(MovingAverages/2):size(FilteredEEGEnvelope,1));
        
        MeanData = mean (X);
        SmoothedFilteredEEGEnvelope (i) = MeanData;
    end

end

clear MovingAverages FilteredEEGEnvelope i X MeanData