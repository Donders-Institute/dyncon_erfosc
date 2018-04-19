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

data_lp = erf_osc_analysis_lcmv_orientation(subj, erfoi, 'reverse', 'lp'); % optimal dipole orientation is chosen based on [0 0.5] after stimulus-reversal.
% data is timelocked to reversal or behavioral response, and low pass
% filtered. 
data = erf_osc_analysis_lcmv_orientation(subj, erfoi, 'reverse', 'no'); % no filtering/project/3011085.02/results/freq/sub-00$i/tzero/gamma_angle_all.mat

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
freq = ft_freqanalysis(cfg, data_lp);

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
inputTime = -inputTime;

%% estimate angle
cfg=[];
cfg.latency = [-1.5+1/fs 0.5];
data_estangle = ft_selectdata(cfg, data);
data_estangle = rmfield(data_estangle, 'cfg');

cfg=[];
cfg.latency = [inputTime inputTime+0.250];
data_erf = ft_selectdata(cfg, data_estangle);

nTrials = length(data_estangle.trial);
nBins = 12;
binSize = floor(nTrials/nBins);

if ischar(freqoistr)
    if strcmp(freqoistr, 'all')
        freqoi = 2:2:100;
    elseif strcmp(freqoistr, 'gamma')
        freqoi = peakFreq_gamma;
    end
end
t = data_estangle.time{1};
m=1;
for frq = freqoi
%
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
%{   
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
fcomp          = ft_freqanalysis(cfg, data_estangle); 
%}
angle_rad   = angle(fcomp.fourierspctrm); % in radians [0, 2pi]
angle_deg   = rad2deg(angle_rad);

%% binning
% estimate angles closest to six center angles. for angles close to 0 and
% 360 degrees, this is a bit more complicated (there is no circular sorting
% option).

trldata = cat(1, data_erf.trial{:});
[cAngles, binAvg{m}, binAngles{m}] = erf_osc_analysis_binangles(trldata, angle_rad, nBins);


%% Cosine fit
% to do: fit a cosine function (statfun_cosinefit)
% compare it with a fit where random angles are joined (avg fit over N(100)
% reps).

cfg = [];
cfg.cosinefit.statistic = 'complex';
design = cAngles;
design = repmat(design, [size(binAvg{m},1), 1]); % repeat for number of time points
[s{m}, ~] = statfun_cosinefit(cfg, binAvg{m}, design);

% Do the same for random distribution
% combine random angles (without replacement; i.e. every trial is used once)
nRep = 100;

for k=1:nRep
    allidx = 1:nTrials;
    for iBin = 1:nBins
        idxtmp = [];
        idxrand = randsample(allidx, binSize);
        randbin(k,iBin,:) = mean(trldata(idxrand,:),1);
        idxtmp = ismember(allidx, idxrand);
        allidx(idxtmp) = [];
    end
end

cfg = [];
cfg.cosinefit.statistic = 'complex';
cfg.parameter = 'trial';
for k=1:nRep
[srand{k}, ~] = statfun_cosinefit(cfg, squeeze(randbin(k,:,:))', design);
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
filename = sprintf('/project/3011085.02/results/freq/sub-%03d/tneg/gamma_angle_%s', subj, freqoistr);

save(fullfile([filename '.mat']), 'angle_rad', 'inputTime','s', 'srand', 'allampstat', 'allanglestat', 'allampstatrand', 'allanglestatrand', 'peakFreq_gamma', 'binAngles');

ft_diary('off')

end


