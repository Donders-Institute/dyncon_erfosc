function erf_osc_analysis_gamma_phase(subj, erfoi, filtdir)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(erfoi)
    erfoi = 'reversal';
end
if nargin<3 && strcmp(erfoi, 'reversal') || isempty(erfoi) && strcmp(erfoi, 'reversal')
    filtdir = 'reverse';
else filtdir=[];
end

%% initiate diary
% ft_diary('on')

%% load data
erf_osc_datainfo;
load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'lcmvData');
load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj), 'peakFreq_gamma');

data = erf_osc_analysis_lcmv_orientation(subj, erfoi, filtdir); % optimal dipole orientation is chosen based on [0 0.5] after stimulus-reversal.
% data is timelocked to reversal or behavioral response, and low pass
% filtered. 

fs=data.fsample;
%% ITC
cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 6;
cfg.foi          = 1:1:24;% analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.1;   % length of time window
cfg.toi          = -1:1/fs:0.6;                  % time window "slides" in steps of 1/fs
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

% compute inter-trial linear coherence (itlc)
itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension

figure;
plot(itc.time, mean(itc.itpc));
title(sprintf('subject %d', subj));


%% estimate gamma angle
cfg=[];
cfg.latency = [-0.25 0.035];
dataShort = ft_selectdata(cfg, data);

cfg            = [];
cfg.taper      = 'hanning';
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 8;
cfg.foilim     = [peakFreq_gamma peakFreq_gamma];
fcomp          = ft_freqanalysis(cfg, dataShort);
gamAngle       = angle(fcomp.fourierspctrm); % in radians

[idx_peak, ~] = find(gamAngle>pi/4 & gamAngle<3*pi/4);
[idx_through, ~] = find(gamAngle>-3*pi/4 & gamAngle<-pi/4);
[idx_up, ~] = find(gamAngle>-pi/4 & gamAngle<pi/4);
[idx_down, ~] = find(gamAngle>3*pi/4 | gamAngle<-3*pi/4);

N = min([length(idx_peak), length(idx_through), length(idx_up), length(idx_down)]);

cfg=[];
cfg.trials = idx_peak(1:N);
cfg.preproc.demean='yes';
cfg.preproc.baselinewindow = [-0.1 0];
peak = ft_selectdata(cfg, data);
cfg.trials = idx_through(1:N);
through = ft_selectdata(cfg, data);
cfg.trials = idx_up(1:N);
up = ft_selectdata(cfg, data);
cfg.trials = idx_down(1:N);
down = ft_selectdata(cfg, data);

figure;
ft_singleplotER([], up, peak, down, through);
xlim([0 0.4])
legend({'up', 'peak', 'down', 'through'})
%% save
% if isPilot
%     filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_angle', subj);
% else
%     filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_angle', subj);
% end
% save(fullfile([filename '.mat']), 'gamAngle');
% ft_diary('off')

end