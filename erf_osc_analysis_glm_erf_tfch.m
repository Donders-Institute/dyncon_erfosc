function erf_osc_analysis_analysis_erf_tfch(subj, isPilot, freqRange, zeropoint, erfoi, doDSS)
% linear regression of peak amplitude over time-frequency (with fixed
% channel) or over frequency-channel (with fixed (avg) time).

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(freqRange)
    freqRange = 'high'; % can be 'high' or 'low'; depending on which frequency range you want to do regression (2-30 Hz or 28-100 Hz)
end
if nargin<4 || isempty(zeropoint)
    zeropoint = 'reversal'; % other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if nargin<5 || isempty(erfoi) % erf of interest
    erfoi = 'reversal'; % other option 'onset': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if nargin<6 || isempty(doDSS)
  doDSS=false;
end

% Initiate Diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/analysis/erf/pilot-%03d/sub-%03d_dss.mat', subj, subj), 'data_dss');
    %load(sprintf('/project/3011085.02/analysis/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    %load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/tfa_%s.mat', subj, zeropoint));
    if strcmp(erfoi, 'reversal')
        if doDSS
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'reversal', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    else
        if doDSS
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'onset', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj));
        end
    end
end
%% if no data_dss cleaning is done beforehand, do this
if ~doDSS
    data = dataClean;

    idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,6)>data.trialinfo(:,5));
    nTrials = length(idxM);
    
    cfg=[];
    cfg.trials = idxM;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    
    % find out which trials have response after end of trial, so you can
    % exclude them
    cfg=[];
    cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
    data_reversal_tmp = ft_redefinetrial(cfg, data);
    
    for iTrial=1:nTrials
        trlLatency(iTrial) = data_reversal_tmp.time{iTrial}(end);
    end
    idx_trials = find(trlLatency'>((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
    idx_trials_invalid = find(trlLatency'<((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
    
    cfg=[];
    cfg.trials = idx_trials;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    clear data_reversal_tmp
    
    if strcmp(erfoi, 'reversal')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
        data = ft_redefinetrial(cfg, data);
    end
    clear dataClean;
else
    data=data_dss;
    clear data_dss;
end

fs=data.fsample;
nTrials = length(data.trial);

%% select p1 amplitude for regression
% find parieto-occipital channel with maximum amplitude in the 60-120 ms
% window.
cfg=[];
cfg.vartrllength=2;
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);
allchans = tlck.label;

t1p1 = nearest(tlck.time, 0.06);
t2p1 = nearest(tlck.time, 0.12);
cfg=[];
cfg.channel = {'MRO', 'MRP', 'MLO', 'MLP', 'MZO', 'MZP'};
cfg.latency = [tlck.time(t1p1) tlck.time(t2p1)];
tlck = ft_selectdata(cfg, tlck);

cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 30;
tlck = ft_preprocessing(cfg, tlck);

cfg=[];
cfg.avgoverrpt = 'yes';
tlck = ft_selectdata(cfg, tlck);

time = tlck.time;
[~, maxchan] = max(abs(mean(tlck.trial,2))); % find channel with max amplitude
% calculate mean over every window to find out which window has the maximum
% mean amplitude (for the maximum channel!). Take the mean amplitude in
% this latency window as regression weight.
maxchanid = tlck.label(maxchan);
maxchanidx = find(strcmp(maxchanid, data.label));
halfwindowlength = 8;

% find window within 60-120ms where amplitude peaks at this channel
i=1;
for ntrl = halfwindowlength+1:length(time)-halfwindowlength
    win(i,:) = [time(ntrl-8), time(ntrl+8)];
    avg(i,1) = mean(tlck.trial(maxchan,ntrl-8:ntrl+8),2);
    i=i+1;
end
[~, window] = max(abs(avg));
lat = win(window,:);

% select parieto-occipital channels in the peak latency
cfg=[];
cfg.latency = lat;
cfg.avgovertime = 'yes';
cfg.channel = {'MRO', 'MRP', 'MLO', 'MLP', 'MZO', 'MZP'};
data = ft_selectdata(cfg, data);

p1amp_all = cat(2,data.trial{:}); % posterior channels only
p1amp_avgtrl = mean(p1amp_all,2); % average p1 amp over trials
multiplier = sign(p1amp_avgtrl); % get the sign of the average
% find a few representative posterior channels that have highest absolute mean amplitude over trials
[~, idx] = sort(multiplier.*p1amp_avgtrl, 'descend'); 
num_p1chans = 4;
p1chans = idx(1:num_p1chans);
p1chans_id = data.label(p1chans);
% representative p1 amplitude is the mean of the *absolute* value of
% representative channels. (Note that absolute here means sign-flipped if
% a channel's average < 0);
p1amp = multiplier(p1chans)'*p1amp_all(p1chans,:)/num_p1chans;

%% Regression p1 amplitude over time-frequency-channel

% load TFA data
load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_tfa_%s_%s.mat', subj, subj, freqRange, zeropoint), 'tfa');
load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_tfa_%s_%s.mat', subj, subj, freqRange, 'onset'), 'baseline');

% select active window of TFA
if strcmp(zeropoint, 'onset')
    cfg=[];
    cfg.latency = [0.125 0.75]; % shortest baseline window is 1 second, shortest reversal 1 second as well
    tfa = ft_selectdata(cfg, tfa);
else
    cfg=[];
    cfg.latency = [-0.75 -0.125];
    tfa = ft_selectdata(cfg, tfa);
    cfg.latency = [-0.9 0.25];
    baseline = ft_selectdata(cfg, baseline);
end
cfg=[];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
tfa = ft_math(cfg, tfa);
baseline = ft_math(cfg, baseline);


nchan = size(tfa.powspctrm, 2);
nfreq = size(tfa.powspctrm, 3);
ntime = size(tfa.powspctrm, 4);
ntrl = size(p1amp,2);

design = [((p1amp-mean(p1amp))./std(p1amp)); ones(1,ntrl)];

cfg=[];
cfg.analysis.statistic = 'beta';

for freq = 1:nfreq
    for ch = 1:nchan
        dat = [squeeze(squeeze(tfa.powspctrm(:,ch,freq,:)))]'; 
        dat = (dat - repmat(mean(dat,2),[1 length(tfa.trialinfo)]))./(repmat(std(dat,[],2),[1 length(tfa.trialinfo)]));
        tmp = statfun_analysis(cfg, dat, design);
        betas_tmp(freq,ch,:,:) = tmp.stat(:,1);
        
        dat_bl = [squeeze(squeeze(baseline.powspctrm(:,ch,freq,:)))]';
        dat_bl = (dat_bl - repmat(mean(dat_bl,2),[1 length(tfa.trialinfo)]))./(repmat(std(dat_bl,[],2),[1 length(tfa.trialinfo)]));
        tmp_bl = statfun_analysis(cfg, dat_bl, design);
        betas_bl_tmp(freq,ch,:,:) = tmp_bl.stat(:,1);        
    end
end

%% Put betas and shuffles in freq-structure
betas = rmfield(tfa, {'powspctrm', 'cfg'});
betas.powspctrm = permute(betas_tmp, [2,1,3]);
betas.dimord = 'chan_freq_time';
betas_bl = rmfield(betas, 'powspctrm');
betas_bl.powspctrm = permute(betas_bl_tmp, [2,1,3]);

%% Save

if isPilot
    filename = sprintf('/project/3011085.02/analysis/erf/pilot-%03d/sub-%03d_analysis_tf_%s_%s_erf_%s', subj, subj, freqRange, zeropoint, erfoi);
else
    filename = sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_analysis_tf_%s_%s_erf_%s', subj, subj, freqRange, zeropoint, erfoi);
end
save(fullfile([filename '.mat']), 'betas', 'betas_bl', 'lat', 'p1amp', 'maxchanid','p1chans_id', '-v7.3');

ft_diary('off')

