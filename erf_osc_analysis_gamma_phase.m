function erf_osc_analysis_gamma_phase(subj, erfoi, freqoistr)
% Estimate the phase of an oscillation with a specific frequency, just
% before the 'input time' (i.e. the time at which ITC is above zero, which
% is an estimate of the time at which the input arrives in the brain). Bin
% the trials according to six center angles and fit a cosine function to
% the average of each bin.
%   input 1: subject number (1 (default) - 33) NOT 10
%   input 2: ERF of interest ('onset' or 'reversal' (default))
%   input 3: frequency of interest (can be 'all' (2:2:100; default), 'gamma' 
%           (gamma peak frequency), or a manually specified number or vector.

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(erfoi)
    erfoi = 'reversal';
end
if nargin<3 || isempty(freqoistr)
    freqoistr = 'all'; % other option is *gamma*, at which the gamma peak frequency is used
end

%% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'lcmvData');
load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj), 'peakFreq_gamma');

data = erf_osc_analysis_lcmv_orientation(subj, erfoi, 'reverse', 'lp'); % optimal dipole orientation is chosen based on [0 0.5] after stimulus-reversal.
% data is timelocked to reversal or behavioral response, and low pass
% filtered. 


fs=data.fsample;
%% ITC
cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 6;
cfg.foi          = 10;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.10;   % length of time window
cfg.toi          = -2:1/fs:1;                  % time window "slides" in steps of 1/fs
freq = ft_freqanalysis(cfg, data);

% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition
itc = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';

F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials

% compute inter-trial phase coherence (itpc) 
itc.itpc      = F./abs(F);         % divide by amplitude  
itc.itpc      = sum(itc.itpc,1);   % sum angles
itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

t2 = nearest(itc.time, 0);
u = mean(itc.itpc(1:t2));
sigma = std(itc.itpc(1:t2));
inputTime = itc.time(find(itc.itpc>u+5*sigma,1));

%% estimate angle
cfg=[];
cfg.latency = [-1.5+1/fs 0.5];
data = ft_selectdata(cfg, data);
data = rmfield(data, 'cfg');


nTrials = length(data.trial);
nBins = 6;
binSize = floor(nTrials/nBins);

% center angles
centerAngle1 = pi/6; % 30 degrees
centerAngle2 = pi/2; % 90 degrees
centerAngle3 = 5*pi/6; % 150 degrees
centerAngle4 = 7*pi/6; % 210 degrees
centerAngle5 = 3*pi/2; % 270 degrees
centerAngle6 = 11*pi/6; % 330 degrees

if ischar(freqoistr)
    if strcmp(freqoistr, 'all')
        freqoi = 2:2:100;
    elseif strcmp(freqoistr, 'gamma')
        freqoi = peakFreq_gamma;
    end
end
t = data.time{1};
m=1;
for frq = freqoi
    %{
    cfg=[];
    cfg.latency = [inputTime-3/frq inputTime-1/fs];
    dataShort = ft_selectdata(cfg, data);
    
    cfg=[];
    cfg.taper = 'hanning';
    cfg.output = 'fourier';
    cfg.method = 'mtmfft';
    cfg.keeptrials = 'yes';
    cfg.foilim = [frq frq];   
    fcomp = ft_freqanalysis(cfg, dataShort);
%}    
cfg            = [];
% use an alpha taper, which concentrates more on the right side of the
% window instead of the center. This way, we can still use mtmconvol, even
% though you the window is centered on the toi (it can not be centered on
% input time because post input time data is needed for that). So from one
% window before the input time, the data near the input time will influence
% the frequency analysis most thanks to the alpha taper.
cfg.taper      = 'alpha'; 
cfg.output     = 'fourier';
cfg.method     = 'mtmconvol';
cfg.keeptrials = 'yes'; 
cfg.keeptapers = 'yes'; % required for mtmconvol i.c.w. fourier output
cfg.foi        = frq; 
cfg.toi        = t(nearest(t,inputTime-0.5/frq-1/fs)); % The window is centered such that
% it includes the data just before input time. The alpha taper makes that
% the estimate is mostly defined by the right side of the window. Use the
% nearest to work around numerical inaccuracies in the time axis (because
% of sample freq)
cfg.t_ftimwin  = 1/frq; 
fcomp          = ft_freqanalysis(cfg, data); 

angle_rad   = angle(fcomp.fourierspctrm)+pi; % in radians [0, 2pi]
angle_deg   = rad2deg(angle_rad);

%% binning
% estimate angles closest to six center angles. for angles close to 0 and
% 360 degrees, this is a bit more complicated (there is no circular sorting
% option).
[val1a, idx1a] = sort(abs(angle_rad-centerAngle1), 'ascend');
[val1b, idx1b] = sort(abs(angle_rad-(centerAngle1 + 2*pi)), 'ascend');
idx1ab = [idx1a; idx1b];
[~, idx1tmp] = sort([val1a; val1b], 'ascend');
idx1 = idx1ab(idx1tmp(1:nTrials));
val1 = angle_rad(idx1);

[~, idx2] = sort(abs(angle_rad-centerAngle2), 'ascend');
val2 = angle_rad(idx2);
[~, idx3] = sort(abs(angle_rad-centerAngle3), 'ascend');
val3 = angle_rad(idx3);
[~, idx4] = sort(abs(angle_rad-centerAngle4), 'ascend');
val4 = angle_rad(idx4);
[~, idx5] = sort(abs(angle_rad-centerAngle5), 'ascend');
val5 = angle_rad(idx5);

[val6a, idx6a] = sort(abs(angle_rad-centerAngle6), 'ascend');
[val6b, idx6b] = sort(abs(angle_rad-(centerAngle6 - 2*pi)), 'ascend');
idx6ab = [idx6a; idx6b];
[~, idx6tmp] = sort([val6a; val6b], 'ascend');
idx6 = idx6ab(idx6tmp(1:nTrials));
val6 = angle_rad(idx6);

binAngles{m} = [val1(1:binSize), val2(1:binSize), val3(1:binSize), val4(1:binSize), val5(1:binSize), val6(1:binSize)];
% estimate avg angle of phase bins
val1tmp=val1;
val1tmp(val1>pi) = val1tmp(val1>pi) + 2*pi; % take care of edges
angle1real = mean(val1tmp(1:binSize));
angle2real = mean(val2(1:binSize));
angle3real = mean(val3(1:binSize));
angle4real = mean(val4(1:binSize));
angle5real = mean(val5(1:binSize));
val6tmp=val6;
val6tmp(val6<pi) = val6tmp(val6<pi) + 2*pi;
angle6real = mean(val6tmp(1:binSize)); % take care of edges

% select trials for each phase bin
trldata = cat(1, data.trial{:});
idx_inputTime = nearest(data.time{1}, inputTime);
bin1 = trldata(idx1(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);
bin2 = trldata(idx2(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);
bin3 = trldata(idx3(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);
bin4 = trldata(idx4(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);
bin5 = trldata(idx5(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);
bin6 = trldata(idx6(1:binSize), idx_inputTime:idx_inputTime+fs*0.25);

alldata_inputTime = trldata(:, idx_inputTime:idx_inputTime+fs*0.25);

%% Cosine fit
% to do: fit a cosine function (statfun_cosinefit)
% compare it with a fit where random angles are joined (avg fit over N(100)
% reps).

data_inputTime_binned = [mean(bin1,1)', mean(bin2,1)', mean(bin3,1)', mean(bin4,1)', mean(bin5,1)', mean(bin6,1)'];

cfg = [];
cfg.cosinefit.statistic = 'complex';
design = [centerAngle1, centerAngle2, centerAngle3, centerAngle4, centerAngle5, centerAngle6];
% design = [angle1real, angle2real, angle3real, angle4real, angle5real, angle6real];
design = design-pi; % statfun_cosinefit can only swallow angular values in the [-pi, pi] regime
design = repmat(design, [size(data_inputTime_binned,1), 1]);
[s{m}, ~] = statfun_cosinefit(cfg, data_inputTime_binned, design);

% Do the same for random distribution
% combine random angles (without replacement; i.e. every trial is used once)
nRep = 100;

for k=1:nRep
allidx = 1:nTrials;

cfg=[];
idx1rand = randsample(allidx, binSize);
rand_inputTime_bin1(k,:) = mean(alldata_inputTime(idx1rand,:),1);

idxtmp = ismember(allidx, idx1rand);
allidx(idxtmp) = [];
idx2rand = randsample(allidx, binSize);
rand_inputTime_bin2(k,:) = mean(alldata_inputTime(idx2rand,:),1);

idxtmp = ismember(allidx, idx2rand);
allidx(idxtmp) = [];
idx3rand = randsample(allidx, binSize);
rand_inputTime_bin3(k,:) = mean(alldata_inputTime(idx3rand,:),1);

idxtmp = ismember(allidx, idx3rand);
allidx(idxtmp) = [];
idx4rand = randsample(allidx, binSize);
rand_inputTime_bin4(k,:) = mean(alldata_inputTime(idx4rand,:),1);

idxtmp = ismember(allidx, idx4rand);
allidx(idxtmp) = [];
idx5rand = randsample(allidx, binSize);
rand_inputTime_bin5(k,:) = mean(alldata_inputTime(idx5rand,:),1);

idxtmp = ismember(allidx, idx5rand);
allidx(idxtmp) = [];
idx6rand = randsample(allidx, binSize);
rand_inputTime_bin6(k,:) = mean(alldata_inputTime(idx6rand,:),1);
end

% rand_inputTime_binned = [rand_inputTime_bin1; rand_inputTime_bin2; rand_inputTime_bin3; rand_inputTime_bin4; rand_inputTime_bin5; rand_inputTime_bin6]';
rand_inputTime_binned = zeros(nRep, size(alldata_inputTime,2), nBins);
cfg = [];
cfg.cosinefit.statistic = 'complex';
cfg.parameter = 'trial';
for k=1:nRep
rand_inputTime_binned(k,:,:) = [rand_inputTime_bin1(k,:); rand_inputTime_bin2(k,:); rand_inputTime_bin3(k,:); rand_inputTime_bin4(k,:); rand_inputTime_bin5(k,:); rand_inputTime_bin6(k,:)]';    
[srand{k}, ~] = statfun_cosinefit(cfg, squeeze(rand_inputTime_binned(k,:,:)), design);
ampstatrand(k,:) = abs(srand{k}.stat);
anglestatrand(k,:) = angle(srand{k}.stat);
end
ampstatrand = mean(ampstatrand,1);
anglestatrand = mean(anglestatrand,1);

allampstat(m,:) = abs(s{m}.stat);
allanglestat(m,:) = angle(s{m}.stat);
allampstatrand(m,:) = ampstatrand;
allanglestatrand(m,:) = anglestatrand;
m=m+1;
end

%% save

if ~ischar(freqoistr)
    if length(freqoistr)>1
        freqoistr = sprintf('%d_%d', freqoistr(1), freqoistr(end));
    else
        freqoistr = num2str(freqoistr);
    end
end
filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_angle_%s', subj, freqoistr);

save(fullfile([filename '.mat']), 'angle_rad', 'inputTime','s', 'srand', 'allampstat', 'allanglestat', 'allampstatrand', 'allanglestatrand', 'peakFreq_gamma', 'binAngles');

ft_diary('off')

end


