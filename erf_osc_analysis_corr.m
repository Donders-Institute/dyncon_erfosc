function erf_osc_analysis_corr(subj, isPilot, correlation, freqRange, zeropoint, erfoi, doDSS, compareQuartile)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(correlation);
    correlation = input('which correlation do you want to compute? choose *gamma_rt*, *gamma_erf*, *gamma_erf_virtualchan* or *amp_tfr*')
end
if nargin<4 || isempty(freqRange)
    freqRange = 'high';
end
if nargin<5 || isempty(zeropoint)
    zeropoint = 'reversal';
end
if nargin<6 || isempty(erfoi)
    erfoi = 'reversal';
end
if nargin<7 || isempty(doDSS)
    doDSS = 0;
end
if nargin<8 || isempty(compareQuartile)
    compareQuartile = false;
end

% initiate diary
ft_diary('on')
erf_osc_datainfo;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamma power - reaction time %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(correlation, 'gamma_rt');
    for subj=allsubs
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj), 'peakFreq_gamma');
        gammaPow{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaPow'); % load gamma power
        gammaPow{subj} = gammaPow{subj}.gammaPow;
        rt{subj} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj)); % load reaction time power
        rt{subj} = rt{subj}.rt;
        
        
        load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
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
        
        trlLatency=[];
        idx_trials = [];
        idx_trials_invalid = [];
        for iTrial=1:nTrials
            trlLatency(iTrial) = data_reversal_tmp.time{iTrial}(end);
        end
        idx_trials = find(trlLatency'>((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
        idx_trials_invalid = find(trlLatency'<((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
        
        cfg=[];
        cfg.trials = idx_trials;
        data = ft_selectdata(cfg, data);
        trialinfo{subj} = data.trialinfo;
        clear data data_reversal_tmp
        jitter{subj} = (trialinfo{subj}(:,5)-trialinfo{subj}(:,4))/1200;
        
        rt{subj} = log(rt{subj});
        rt{subj} = rt{subj}-mean(rt{subj});
        jitter{subj} = log(jitter{subj});
        jitter{subj} = jitter{subj}-mean(jitter{subj});
        [r_gamma_rt(subj) p_gamma_rt(subj)] = partialcorr(gammaPow{subj}', rt{subj},jitter{subj}, 'type', 'spearman');
        [r_jitter_gamma(subj) p_jitter_gamma(subj)] = partialcorr(gammaPow{subj}', jitter{subj},rt{subj}, 'type', 'spearman');
        [r_jitter_rt(subj) p_jitter_rt(subj)] = partialcorr(rt{subj}, jitter{subj},gammaPow{subj}', 'type', 'spearman');
        
    end
    
    % statistics
    stat_gamma_rt=[];
    [stat_gamma_rt.h stat_gamma_rt.p]= ttest(r_gamma_rt);
    stat_jitter_gamma=[];
    [stat_jitter_gamma.h stat_jitter_gamma.p] = ttest(r_jitter_gamma);
    stat_jitter_rt=[];
    [stat_jitter_rt.h stat_jitter_rt.p] = ttest(r_jitter_rt);
    
elseif strcmp(correlation, 'amp_tfr') || strcmp(correlation, 'gamma_erf')
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
        if strcmp(correlation, 'amp_tfr')
            load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s_%s.mat', subj, freqRange, zeropoint), 'tfa');
            load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s_%s.mat', subj, freqRange, 'onset'), 'baseline');
        elseif strcmp(correlation, 'gamma_erf')
            load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaPow');
        end
    end
    %% if no data_dss cleaning is done beforehand, do this
    if ~doDSS
        data=dataClean;
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
        elseif strcmp(erfoi, 'motor')
            cfg=[];
            cfg.offset = -(data.trialinfo(:,6)-data.trialinfo(:,4));
            data=ft_redefinetrial(cfg, data);
        end
        clear dataClean data_reversal_tmp trlLatency;
    else
        data=data_dss;
        clear data_dss;
    end
    
    fs=data.fsample;
    nTrials = length(data.trial);
    
    
    %% compute P1 amplitude
     cfg=[];
        cfg.lpfilter = 'yes';
        cfg.lpfilttype = 'firws';
        cfg.lpfreq = 30;
        data = ft_preprocessing(cfg, data);
    if strcmp(correlation, 'amp_tfr')
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
        for nTrials = halfwindowlength+1:length(time)-halfwindowlength
            win(i,:) = [time(nTrials-8), time(nTrials+8)];
            avg(i,1) = mean(tlck.trial(maxchan,nTrials-8:nTrials+8),2);
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
        
        nTrials = length(p1amp);
        
        
        %% preprocess TFR data
        % select active window of TFR
        if strcmp(zeropoint, 'onset')
            cfg=[];
            cfg.latency = [0.125 0.75]; % shortest baseline window is 1 second, shortest reversal 1 second as well
            tfa = ft_selectdata(cfg, tfa);
        else
            cfg=[];
            cfg.latency = [-0.5 -0.125];
            tfa = ft_selectdata(cfg, tfa);
            cfg.latency = [-0.650 0.25];
            baseline = ft_selectdata(cfg, baseline);
        end
        
        
        nchan = size(tfa.powspctrm, 2);
        nfreq = size(tfa.powspctrm, 3);
        ntime = size(tfa.powspctrm, 4);
        
        
        if compareQuartile
            nQ = 4; % number of quartiles
            qSize = floor(nTrials/4); % quartile size
            [val idx] = sort(p1amp, 2, 'descend');
            q1 = idx(1:qSize); % trials with highest p1 amp
            q4 = idx(end-(qSize-1):end); % trials with lowest p1 amp
            
            cfg=[];
            cfg.avgoverrpt = 'yes';
            cfg.trials = q1;
            tfa_q1 = ft_selectdata(cfg, tfa);
            cfg.trials = q4;
            tfa_q4 = ft_selectdata(cfg, tfa);
            
            cfg=[];
            cfg.parameter = 'powspctrm';
            cfg.operation = 'log10';
            tfa_q1 = ft_math(cfg, tfa_q1);
            tfa_q4 = ft_math(cfg, tfa_q4);
        else
            cfg=[];
            cfg.parameter = 'powspctrm';
            cfg.operation = 'log10';
            tfa = ft_math(cfg, tfa);
            baseline = ft_math(cfg, baseline);
            
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
    elseif strcmp(correlation, 'gamma_erf')
        nTrials = length(data.trial);
        
        if compareQuartile
            nQ = 4; % number of quartiles
            qSize = floor(nTrials/4); % quartile size
            [val idx] = sort(gammaPow, 2, 'descend');
            q1 = idx(1:qSize); % trials with highest p1 amp
            q4 = idx(end-(qSize-1):end); % trials with lowest p1 amp
                      
            cfg=[];
            cfg2=[];
            cfg.trials = q1;
            data_q1 = ft_selectdata(cfg, data);
            cfg2.vartrllength = 2;
            tlck_q1 = ft_timelockanalysis(cfg2, data_q1);
            cfg.trials = q4;
            data_q4 = ft_selectdata(cfg, data);
            tlck_q4 = ft_timelockanalysis(cfg2, data_q4);
        end
        
    end
elseif strcmp(correlation, 'gamma_erf_virtualchan')
        load(sprintf('/project/3011085.02/results/freq/sub-%03d/allori/gamma_virtual_channel.mat', subj), 'gammaPow');
        lcmvData = erf_osc_analysis_lcmv_orientation(subj, erfoi); 

end
    %% save
    if strcmp(correlation, 'gamma_rt');
        filename = '/project/3011085.02/results/stat_partcorr_gamma_rt';
        save(fullfile([filename '.mat']), 'stat_gamma_rt', 'r_gamma_rt', 'p_gamma_rt', 'rt', 'gammaPow', 'r_jitter_gamma', 'p_jitter_gamma', 'stat_jitter_gamma', 'r_jitter_rt', 'p_jitter_rt', 'stat_jitter_rt');
    elseif strcmp(correlation, 'amp_tfr') && ~compareQuartile
        filename = sprintf('/project/3011085.02/results/freq/sub-%03d/corr_amp_tfr', subj);
        save(fullfile([filename '.mat']), 'z_act', 'z_bl', 'r_act', 'r_bl', 'lat', 'p1amp', 'maxchanid','p1chans_id');
    elseif strcmp(correlation, 'amp_tfr') && compareQuartile
        filename = sprintf('/project/3011085.02/results/freq/sub-%03d/corr_amp_tfr_quartile', subj);
        save(fullfile([filename '.mat']), 'tfa_q1', 'tfa_q4', 'lat', 'p1amp', 'maxchanid','p1chans_id', 'val', 'idx', 'q1', 'q4');
    elseif strcmp(correlation, 'gamma_erf') && compareQuartile
        filename = sprintf('/project/3011085.02/results/erf/sub-%03d/corr_gamma_erf_quartile_%s', subj, erfoi);
        save(fullfile([filename '.mat']), 'tlck_q1', 'tlck_q4', 'gammaPow', 'val', 'idx', 'q1', 'q4');
    end
    ft_diary('off')

