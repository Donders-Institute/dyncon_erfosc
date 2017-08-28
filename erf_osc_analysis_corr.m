function erf_osc_analysis_corr(subj, isPilot, correlation, freqRange, zeropoint, erfoi, doDSS)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(correlation);
    correlation = input('which correlation do you want to compute? choose *gamma_rt* or *amp_tfr*')
end
if nargin<4 || isempty(freqRange)
    freqRange = 'high';
end
if nargin<5 || isempty(zeropoint)
    zeropoint = 'reversal';
end
if nargin<6 || isempty(erfoi)
    erfoi = 'rerversal';
end
if nargin<7 || isempty(doDSS)
    doDSS = 0;
end

% initiate diary
ft_diary('on')
erf_osc_datainfo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamma power - reaction time %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(correlation, 'gamma_rt');
    for subj=allsubs
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_peak.mat', subj), 'peakFreq_gamma');
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan'); % load gamma power
        rt{subj} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj)); % load gamma power
        rt{subj} = rt{subj}.rt;
        
        for i=1:length(gammaChan.trial)
            gammaPow_tmp(i) = log(gammaChan.trial(i).pow);
        end
        gammaPow{subj} = gammaPow_tmp-mean(gammaPow_tmp);
        clear gammaChan gammaPow_tmp
        [r(subj) p(subj)] = corr(gammaPow{subj}', rt{subj}, 'type', 'spearman');
    end 
    
    % statistics
    stat=[];
    [stat.h stat.p]= ttest(r);
    
elseif strcmp(correlation, 'amp_tfr')
    %%%%%%%%%%%%%%%%%%%%%%
    % P1 amplitude - TFR %
    %%%%%%%%%%%%%%%%%%%%%%
    %% load data
    if isPilot
        load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    else
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
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s_%s.mat', subj, freqRange, zeropoint), 'tfa');
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s_%s.mat', subj, freqRange, 'onset'), 'baseline');
    end
    %% if no data_dss cleaning is done beforehand, do this
    if ~doDSS
        data=dataClean;
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
        clear dataClean;
    else
        data=data_dss;
        clear data_dss;
    end
    
    fs=data.fsample;
    nTrials = length(data.trial);
    
    %% compute P1 amplitude
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
    i=1;
    for ntrl = halfwindowlength+1:length(time)-halfwindowlength
        win(i,:) = [time(ntrl-8), time(ntrl+8)];
        avg(i,1) = mean(tlck.trial(maxchan,ntrl-8:ntrl+8),2);
        i=i+1;
    end
    [~, window] = max(abs(avg));
    lat = win(window,:);
    
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
    % standardize p1 data
    p1amp = (p1amp-mean(p1amp,2))./std(p1amp,[],2);
    
    %% preprocess TFR data
    % select active window of TFR
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
    
    %% Correlation
    for freq = 1:nfreq
        for ch = 1:nchan
            [r_act_tmp, ~] = corr(p1amp',  squeeze(squeeze(tfa.powspctrm(:,ch,freq,:))), 'type', 'spearman');
            [r_bl_tmp, ~] = corr(p1amp',  squeeze(squeeze(baseline.powspctrm(:,ch,freq,:))), 'type', 'spearman');
            r_act(ch,freq,:) = r_act_tmp;
            r_bl(ch,freq,:) = r_bl_tmp;
        end
    end
    
    % Fisher z-transformation
    z_act_tmp = atanh(r_act);
    z_bl_tmp = atanh(r_bl);
    
    %% put correlations in FieldTrip freq structure
    z_act = rmfield(tfa, {'powspctrm', 'cfg'});
    z_act.dimord = 'chan_freq_time';
    z_bl = z_act;
    z_act.powspctrm = z_act_tmp;
    z_bl.powspctrm = z_bl_tmp;
    
    
end
%% save
if strcmp(correlation, 'gamma_rt');
    filename = '/project/3011085.02/results/stat_corr_gamma_rt';
    save(fullfile([filename '.mat']), 'stat', 'r', 'p', 'rt', 'gammaPow');
elseif strcmp(correlation, 'amp_tfr')
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/corr_amp_tfr', subj);
    save(fullfile([filename '.mat']), 'z_act', 'z_bl', 'r_act', 'r_bl', 'lat', 'p1amp', 'maxchanid','p1chans_id');
end
ft_diary('off')
end
