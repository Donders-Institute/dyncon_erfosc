function erf_osc_analysis_glm_erf_tfch(subj, isPilot, freqRange, zeropoint, erfoi, doDSS)
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
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'reversal', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    else
        if doDSS
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'onset', false);
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
else
    data=data_dss;
end

fs=data.fsample;
nTrials = length(data.trial);


%% select p1 window for regression
cfg=[];
cfg.vartrllength=2;
tlck = ft_timelockanalysis(cfg, data);

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
maxchanidx = find(strcmp(maxchanid, data.label));
halfwindowlength = 8;
i=1;
for ntrl = halfwindowlength+1:length(time)-halfwindowlength;
    win(i,:) = [time(ntrl-8), time(ntrl+8)];
    avg(i,1) = mean(tlck.avg(maxchan,ntrl-8:ntrl+8),2);
    i=i+1;
end
[~, window] = max(abs(avg));
lat = win(window,:);


%% Regression p1 amplitude over time-frequency-channel
cfg=[];
cfg.latency = [-1.5 0.65];
data = ft_selectdata(cfg, data);
trialdata = cat(3,data.trial{:});

idxtime = [nearest(data.time{1}, lat(1)) : nearest(data.time{1}, lat(2))];
p1amp = squeeze(mean(trialdata(:,idxtime,:),2));

% baselinecorrect with average baseline over trials
if strcmp(freqRange, 'high');
    if strcmp(zeropoint, 'onset')
        cfg=[];
        cfg.latency = [-1 1]; % shortest baseline window is 1 second, shortest reversal 1 second as well
        tfaHigh = ft_selectdata(cfg, tfaHigh);
    else
        cfg=[];
        cfg.latency = [-1.5 0];
        tfaHigh = ft_selectdata(cfg, tfaHigh);
    end
    
    for freq=1:19
        for ch=1:length(data.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_h = squeeze(squeeze(tfaHigh.powspctrm(:,ch,freq,:)));
            betas(freq,ch,:,:) = design'\Y_h;
        end
    end
    
    nchan = size(tfaHigh.powspctrm, 2);
    nfreq = size(tfaHigh.powspctrm, 3);
    ntime = size(tfaHigh.powspctrm, 4);
    ntrl = size(p1amp,2);
    design_bias = ones(1, ntrl);
    p1amp_norm = (p1amp - repmat(mean(p1amp,2), [1 ntrl]))./repmat(std(p1amp, [], 2), [1, ntrl]);
    
    numShuffles=1000;
    avg_shuffles = zeros(nchan, nfreq, ntime);
    std_shuffles = zeros(nchan, nfreq, ntime);
    for freq=1:19
        Y1 = squeeze(tfaHigh.powspctrm(:,:,freq,:));
        for ch=1:nchan;
            Y2 = squeeze(Y1(:,ch,:));
            design = [design_bias; p1amp_norm(ch,:)];
            betas_shuffles=zeros(numShuffles, ntime);
            for iShuffle=1:numShuffles
                design2 = design(:,randperm(ntrl));
                tmp = design2'\Y2;
                betas_shuffles(iShuffle, :) = tmp(2,:);
            end
            avg_shuffles(ch, freq, :) = mean(betas_shuffles,1);
            std_shuffles(ch, freq, :) = std(betas_shuffles,[],1);
        end
    end
    
elseif strcmp(freqRange, 'low')
    if strcmp(zeropoint, 'onset')
        cfg=[];
        cfg.latency = [-1 1]; % shortest baseline window is 1 second
        tfaLow = ft_selectdata(cfg, tfaLow);
    else
        cfg=[];
        cfg.latency = [-1.5 0];
        tfaLow = ft_selectdata(cfg, tfaLow);
    end
    
    for freq=1:15
        for ch=1:length(data.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_l = squeeze(squeeze(tfaLow.powspctrm(:,ch,freq,:)));
            betas(freq,ch,:,:) = design'\Y_l;
        end
    end
    
    nchan = size(tfaLow.powspctrm, 2);
    nfreq = size(tfaLow.powspctrm, 3);
    ntime = size(tfaLow.powspctrm, 4);
    ntrl = size(p1amp,2);
    design_bias = ones(1, ntrl);
    p1amp_norm = (p1amp - repmat(mean(p1amp,2), [1 ntrl]))./repmat(std(p1amp, [], 2), [1, ntrl]);
    
    numShuffles=1000;
    avg_shuffles = zeros(nfreq, nchan, ntime);
    std_shuffles = zeros(nfreq, nchan, ntime);
    for freq=1:nfreq
        Y1 = squeeze(tfaLow.powspctrm(:,:,freq,:));
        for ch=1:nchan;
            Y2 = squeeze(Y1(:,ch,:));
            design = [design_bias; p1amp_norm(ch,:)];
            betas_shuffles=zeros(numShuffles, ntime);
            for iShuffle=1:numShuffles
                design2 = design(:,randperm(ntrl));
                tmp = design2'\Y2;
                betas_shuffles(iShuffle, :) = tmp(2,:);
            end
            avg_shuffles(ch, freq, :) = mean(betas_shuffles,1);
            std_shuffles(ch, freq, :) = std(betas_shuffles,[],1);
        end
    end
end


%% planar gradiant transformation of beta weights
% make timelocked structure where planar gradient transformation can be
% applied to (that's why dimord is strange)
if strcmp(freqRange, 'high')
    tl=[];
    tl.avg    = squeeze(betas(:,:,2,:));
    tl.time   = tfaHigh.time;
    tl.dimord = 'subj_chan_time';
    tl.label  = tfaHigh.label;
    tl.grad   = tfaHigh.grad;
    
    % also put betas in a time-freq structure
    tf = rmfield(tl,'avg');
    tf.dimord = 'chan_freq_time';
    tf.powspctrm = permute(tl.avg, [2,1,3]);
    tf.freq = tfaHigh.freq;
    
    % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
    cfg.planarmethod    = 'sincos';
    tlPlanar           = ft_megplanar(cfg, tl);
    
    cfg           = [];
    bPlanarCmb  = ft_combineplanar(cfg,tlPlanar);
    
    % put it back in a freq-data structure
    bPlanarCmb.powspctrm = permute(bPlanarCmb.trial, [2,1,3]);
    bPlanarCmb           = rmfield(bPlanarCmb, 'trial');
    bPlanarCmb.freq      = tfaHigh.freq;
    bPlanarCmb.dimord    = 'chan_freq_time';
    
    tl_shuffle_avg      = rmfield(tl, 'avg');
    tl_shuffle_avg.avg  = permute(avg_shuffles, [2,1,3]);
    tl_shuffle_std      = rmfield(tl, 'avg');
    tl_shuffle_std.avg  = permute(std_shuffles,[2,1,3]);
    
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tl_shuffle_avg);
    cfg.planarmethod    = 'sincos';
    tl_shuffle_avgPlanar= ft_megplanar(cfg, tl_shuffle_avg);
    tl_shuffle_stdPlanar= ft_megplanar(cfg, tl_shuffle_std);
    
    cfg                      = [];
    tl_shuffle_avgCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_avgPlanar);
    tl_shuffle_stdCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_stdPlanar);
    
    shuffle_avgCmbPl = rmfield(bPlanarCmb, 'powspctrm');
    shuffle_avgCmbPl.powspctrm = permute(tl_shuffle_avgCmbPlanar.trial, [2,1,3]);
    shuffle_stdCmbPl = rmfield(bPlanarCmb, 'powspctrm');
    shuffle_stdCmbPl.powspctrm = permute(tl_shuffle_stdCmbPlanar.trial, [2,1,3]);
    
elseif strcmp(freqRange, 'low')
    % Do the same for low frequencies.
    tl=[];
    tl.avg    = squeeze(betas(:,:,2,:));
    tl.time   = tfaLow.time;
    tl.dimord = 'subj_chan_time';
    tl.label  = tfaLow.label;
    tl.grad   = tfaLow.grad;
    
    % also put betas in a time-freq structure
    tf = rmfield(tl,'avg');
    tf.dimord = 'chan_freq_time';
    tf.powspctrm = permute(tl.avg, [2,1,3]);
    tf.freq = tfaLow.freq;
    
    % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
    cfg.planarmethod    = 'sincos';
    tlPlanar           = ft_megplanar(cfg, tl);
    
    cfg           = [];
    bPlanarCmb  = ft_combineplanar(cfg,tlPlanar);
    
    % put it back in a freq-data structure
    bPlanarCmb.powspctrm = permute(bPlanarCmb.trial, [2,1,3]);
    bPlanarCmb           = rmfield(bPlanarCmb, 'trial');
    bPlanarCmb.freq      = tfaLow.freq;
    bPlanarCmb.dimord    = 'chan_freq_time';
    
    tl_shuffle_avg      = rmfield(tl, 'avg');
    tl_shuffle_avg.avg  = permute(avg_shuffles, [2,1,3]);
    tl_shuffle_std      = rmfield(tl, 'avg');
    tl_shuffle_std.avg  = permute(std_shuffles,[2,1,3]);
    
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tl_shuffle_avg);
    cfg.planarmethod    = 'sincos';
    tl_shuffle_avgPlanar= ft_megplanar(cfg, tl_shuffle_avg);
    tl_shuffle_stdPlanar= ft_megplanar(cfg, tl_shuffle_std);
    
    cfg                      = [];
    tl_shuffle_avgCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_avgPlanar);
    tl_shuffle_stdCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_stdPlanar);
    
    shuffle_avgCmbPl = rmfield(bPlanarCmb, 'powspctrm');
    shuffle_avgCmbPl.powspctrm = permute(tl_shuffle_avgCmbPlanar.trial, [2,1,3]);
    shuffle_stdCmbPl = rmfield(bPlanarCmb, 'powspctrm');
    shuffle_stdCmbPl.powspctrm = permute(tl_shuffle_stdCmbPlanar.trial, [2,1,3]);
end

%% Save

if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_tf_%s_%s_erf_%s', subj, freqRange, zeropoint, erfoi);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_%s_erf_%s', subj, freqRange, zeropoint, erfoi);
end

save(fullfile([filename '.mat']), 'betas','bPlanarCmb','tf','lat','maxchanid','maxchanidx', 'shuffle_avgCmbPl', 'shuffle_stdCmbPl', 'p1amp', '-v7.3');



ft_diary('off')

