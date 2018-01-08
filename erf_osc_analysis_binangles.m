function [centerAngles, binAvg, binAngles] = erf_osc_analysis_binangles(data, angles, nBins)


[nTrials, nTime] = size(data); % data should be trial x time
binSize = floor(nTrials/nBins);

centerAngles = -pi:2*pi/(nBins):pi-2*pi/(nBins);
binAvg = zeros(nBins, nTime);
for iBin = 1:nBins
    % turn the center angle to zero
    fdata = exp(1i*angles).*exp(-1i*centerAngles(iBin));
    angle_tozero = angle(fdata);
    
    [~, idx] = sort(abs(angle_tozero), 'ascend');
    val = angles(idx(1:binSize));
    idx = idx(1:binSize);
    
    binAngles(:, iBin) = val;
    binData = data(idx,:);
    binAvg(iBin,:) = mean(binData,1);
end
    binAvg = binAvg';
    



