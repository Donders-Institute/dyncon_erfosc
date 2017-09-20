function erf_osc_analysis_glm_gamma_time(subj, isPilot, erfoi, doDSS)
% do a linear regression of pre-change gamma power over time.
if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(erfoi)
    erfoi = 'reversal'; % can be *onset*, *reversal*, *motor*
end
if nargin<4 || isempty(doDSS)
    doDSS = false;
end


% Initiate Diary
ft_diary('on')


%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
    if doDSS
        [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'reversal', false);
    else
        load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        data = dataClean;
        clear dataClean
    end
end
fs=data.fsample;
if ~doDSS
    
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
    
    
    if strcmp(erfoi, 'reversal')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
        data=ft_redefinetrial(cfg, data);
    elseif strcmp(erfoi, 'motor')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,6)-data.trialinfo(:,4));
        data=ft_redefinetrial(cfg, data);
    end
    clear data_reversal_tmp trlLatency
end

for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
end
gammaPow = (gammaPow-mean(gammaPow))/std(gammaPow);
nTrials = length(data.trial);

[~, idxMax] = sort(gammaPow, 2, 'descend');

%% GLM on all trials
% first select data epochs. Filtering the data might introduce past data to
% future when using a forward filter). Therefore apply reverse filter.
% filter before data-cutting to avoid edge-effects.

cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 30;
cfg.lpfiltdir = 'onepass-reverse-zerophase';
data = ft_preprocessing(cfg, data);

% now cut out the segment of interest.
cfg=[];
if strcmp(erfoi, 'motor')
    cfg.latency = [-0.5 0];
    active = ft_selectdata(cfg, data);
else
    cfg.latency = [0 0.5];
    active = ft_selectdata(cfg, data);
    cfg.latency = [-0.5 0];
    baseline = ft_selectdata(cfg, data);
end

active.trial = cat(3,active.trial{:});
active.trial = permute(active.trial, [3,1,2]);
active.time = active.time{1};

if ~strcmp(erfoi, 'motor')
    baseline.trial = cat(3,baseline.trial{:});
    baseline.trial = permute(baseline.trial, [3,1,2]);
    baseline.time = active.time;
end

design = [gammaPow; ones(size(gammaPow))];

cfg=[];
cfg.glm.statistic = 'beta';

for k=1:length(active.label)
    dat = [squeeze(active.trial(:,k,:))]';
    dat = (dat - repmat(mean(dat,2),[1 length(data.trialinfo)]))./(repmat(std(dat,[],2),[1 length(data.trialinfo)]));
    tmp = statfun_glm(cfg, dat, design);
    betas_tmp(k,:) = tmp.stat(:,1);
    
    if ~strcmp(erfoi, 'motor')
        dat_bl = [squeeze(baseline.trial(:,k,:))]';
        dat_bl = (dat_bl - repmat(mean(dat_bl,2),[1 length(data.trialinfo)]))./(repmat(std(dat_bl,[],2),[1 length(data.trialinfo)]));
        tmp_bl = statfun_glm(cfg, dat_bl, design);
        betas_bl_tmp(k,:) = tmp_bl.stat(:,1);
    end
end

% put beta weights in timelock structure
betas        = rmfield(data,{'trial', 'cfg'});
betas.avg    = betas_tmp;
betas.time   = active.time;
betas.dimord = 'chan_time';
if ~strcmp(erfoi, 'motor')
    betas_bl     = rmfield(betas, 'avg');
    betas_bl.avg = betas_bl_tmp;
end

cfg                 = [];
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, betas);
cfg.planarmethod    = 'sincos';
betas_planar        = ft_megplanar(cfg, betas);
betas_plcmb         = ft_combineplanar([], betas_planar);
if ~strcmp(erfoi, 'motor')
    betas_bl_planar     = ft_megplanar(cfg, betas_bl);
    betas_bl_plcmb      = ft_combineplanar([], betas_bl_planar);
else
    betas_bl_plcmb = 'use baseline in erfoi = reversal';
    betas_bl = 'use baseline in erfoi = reversal';
end

%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_gamma_time_%s', subj, erfoi);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time_%s', subj, erfoi);
end
save(fullfile([filename '.mat']), 'betas_plcmb','betas_bl_plcmb', 'betas','betas_bl', '-v7.3');
ft_diary('off')

