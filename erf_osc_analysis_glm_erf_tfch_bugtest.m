function erf_osc_analysis_glm_erf_tfch_bugtest(subj, isPilot, freqRange, zeropoint, erfoi, doDSS)
% linear regression of peak amplitude over time-frequency (with fixed
% channel) or over frequency-channel (with fixed (avg) time).

if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = false;
end
if isempty(isPilot);
    isPilot = false;
end
if nargin<3
    freqRange = 'high'; % can be 'high' or 'low'; depending on which frequency range you want to do regression (2-30 Hz or 28-100 Hz)
end
if isempty(freqRange)
    freqRange = 'high';
end
if nargin<4
    zeropoint = 'onset'; % other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if isempty(zeropoint);
    zeropoint = 'onset';
end
if nargin<5 % erf of interest
    erfoi = 'reversal'; % other option 'onset': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if isempty(erfoi);
    erfoi = 'reversal';
end

% Initiate Diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s.mat', subj, zeropoint));
    if strcmp(erfoi, 'reversal')
        if doDSS
            [data_dss, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'reversal', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    else
        if doDSS
            [data_dss, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'onset', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    end
end
%% if no data_dss cleaning is done beforehand, do this
if ~doDSS
    data=dataClean;
    cfg=[];
    idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
    nTrials = length(idxM);
    
    cfg=[];
    cfg.trials = idxM;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    
    if strcmp(erfoi, 'reversal')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
        data = ft_redefinetrial(cfg, data);
    end
    data_dss=data;
end

fs=data_dss.fsample;
nTrials = length(data_dss.trial);

%% select p1 window for regression
cfg=[];
cfg.vartrllength=2;
tlck = ft_timelockanalysis(cfg, data_dss);

t1p1 = nearest(tlck.time, 0.06);
t2p1 = nearest(tlck.time, 0.12);
cfg=[];
cfg.channel = {'MRO', 'MRP', 'MLO', 'MRO', 'MZO', 'MZP'};
cfg.latency = [tlck.time(t1p1) tlck.time(t2p1)];
tlck = ft_selectdata(cfg, tlck);

time = tlck.time;
[~, maxchan] = max(abs(mean(tlck.avg,2))); % find channel with max amplitude
% calculate mean over every window to find out which window has the maximum
% mean amplitude (for the maximum channel!). Take the mean amplitude in
% this latency window as regression weight.
maxchanid = tlck.label(maxchan);
halfwindowlength = 8;
i=1;
for t = halfwindowlength+1:length(time)-halfwindowlength;
    win(i,:) = [time(t-8), time(t+8)];
    avg(i,1) = mean(tlck.avg(maxchan,t-8:t+8),2);
    i=i+1;
end
[~, window] = max(abs(avg));
lat = win(window,:);

%% Regression p1 amplitude over time-frequency-channel
cfg=[];
cfg.latency = [-1.5 0.65];
data_dss = ft_selectdata(cfg, data_dss);
trialdata = cat(3,data_dss.trial{:});

idxtime = [nearest(data_dss.time{1}, lat(1)) : nearest(data_dss.time{1}, lat(2))];
p1amp = squeeze(mean(trialdata(:,idxtime,:),2));

% baselinecorrect with average baseline over trials
if strcmp(freqRange, 'high');
    if strcmp(zeropoint, 'onset')
        cfg=[];
        cfg.latency = [-1 1.75]; % shortest baseline window is 1 second
        tfaHigh = ft_selectdata(cfg, tfaHigh);
    else
        cfg=[];
        cfg.latency = [-1.5 0.5];
        tfaHigh = ft_selectdata(cfg, tfaHigh);
    end
     load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_onset.mat', subj), 'baselineH');
      baselineH.time = tfaHigh.time;
      baselineH.dimord = tfaHigh.dimord;
      baselineH.powspctrm = repmat(baselineH.powspctrm, [1,1,length(baselineH.time), size(tfaHigh.powspctrm, 1)]);
    for freq=1:19
        for ch=1:length(data_dss.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_h = squeeze(squeeze(tfaHigh.powspctrm(:,ch,freq,:)));
            betas_h(freq,ch,:,:) = design'\Y_h;
        end
    end
    
elseif strcmp(freqRange, 'low')
    if strcmp(zeropoint, 'onset')
        cfg=[];
        cfg.latency = [-1 1.75]; % shortest baseline window is 1 second
        tfaLow = ft_selectdata(cfg, tfaLow);
    end
    
    for freq=1:15
        for ch=1:length(data_dss.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_l = squeeze(squeeze(tfaLow.powspctrm(:,ch,freq,:)));
            betas_l(freq,ch,:,:) = design'\Y_l;
        end
    end
end

%% planar gradiant transformation of beta weights
% make timelocked structure where planar gradient transformation can be
% applied to (that's why dimord is strange)
if strcmp(freqRange, 'high')
    tlh=[];
    tlh.avg    = squeeze(betas_h(:,:,2,:));
    tlh.time   = tfaHigh.time;
    tlh.dimord = 'subj_chan_time';
    tlh.label  = tfaHigh.label;
    tlh.grad   = tfaHigh.grad;
    
    % also put betas in a time-freq structure
    tfh = rmfield(tlh,'avg');
    tfh.dimord = 'chan_freq_time';
    tfh.powspctrm = permute(tlh.avg, [2,1,3]);
    tfh.freq = tfaHigh.freq;
    
    % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tlh);
    cfg.planarmethod    = 'sincos';
    tlhPlanar           = ft_megplanar(cfg, tlh);
    
    cfg           = [];
    bhPlanarCmb  = ft_combineplanar(cfg,tlhPlanar);
    
    % put it back in a freq-data structure
    bhPlanarCmb.powspctrm = permute(bhPlanarCmb.trial, [2,1,3]);
    bhPlanarCmb           = rmfield(bhPlanarCmb, 'trial');
    bhPlanarCmb.freq      = tfaHigh.freq;
    bhPlanarCmb.dimord    = 'chan_freq_time';
    
elseif strcmp(freqRange, 'low')
    % Do the same for low frequencies.
    tll=[];
    tll.avg    = betas_l;
    tll.time   = tfaLow.time;
    tll.dimord = 'subj_chan_time';
    tll.label  = tfaLow.label;
    tll.grad   = tfaLow.grad;
    
    % also put betas in a time-freq structure
    tfl = rmfield(tll,'avg');
    tfl.dimord = 'chan_freq_time';
    tfl.powspctrm = permute(tll.avg, [2,1,3]);
    tfl.freq = tfaLow.freq;
    
    % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tll);
    cfg.planarmethod    = 'sincos';
    tllPlanar           = ft_megplanar(cfg, tll);
    
    cfg           = [];
    blPlanarCmb  = ft_combineplanar(cfg,tllPlanar);
    
    % put it back in a freq-data structure
    blPlanarCmb.powspctrm = permute(blPlanarCmb.trial, [2,1,3]);
    blPlanarCmb           = rmfield(blPlanarCmb, 'trial');
    blPlanarCmb.freq      = tfaLow.freq;
    blPlanarCmb.dimord    = 'chan_freq_time';
end

%% Save

if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_tf_%s_%s_erf_%s_bugtest', subj, freqRange, zeropoint, erfoi);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_%s_erf_%s_bugtest', subj, freqRange, zeropoint, erfoi);
end
if strcmp(freqRange, 'high')
    save(fullfile([filename '.mat']), 'betas_h','bhPlanarCmb','tfh','lat','maxchanid', '-v7.3');
elseif strcmp(freqRange, 'low')
    save(fullfile([filename '.mat']), 'betas_l','bhPlanarCmb','tfl', 'lat','maxchanid', '-v7.3');
end

ft_diary('off')

