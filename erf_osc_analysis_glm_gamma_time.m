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
    elseif strcmp(erfoi, 'motor')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,6)-data.ltrialinfo(:,4));
        data=ft_redefinetrial(cfg, data);
    end
end

for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
end
gammaPow = (gammaPow-mean(gammaPow))/std(gammaPow);
nTrials = length(data.trial);

[~, idxMax] = sort(gammaPow, 2, 'descend');

%% GLM on all trials

cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 30;
% cfg.lpfiltdir = 'onepass-reverse';
data = ft_preprocessing(cfg, data);
if ~doDSS
    cfg=[];
    cfg.latency = [-1 0.65];
    data = ft_selectdata(cfg, data);
end

data_conc = data;
data_conc.trial = cat(3,data_conc.trial{:});
data_conc.trial = permute(data_conc.trial, [3,1,2]);
data_conc.time = data_conc.time{1};
if ~strcmp(erfoi,'motor')
    cfg=[];
    cfg.latency = [-0.5 0];
    baseline = ft_selectdata(cfg, data_conc);
end
cfg.latency = [0 0.5];
active = ft_selectdata(cfg, data_conc);


% design = [gammaPow zeros(size(gammaPow)); zeros(size(gammaPow)) gammaPow; ones(size(gammaPow)) ones(size(gammaPow)); ...
%     0.5*ones(size(gammaPow)) -0.5*ones(size(gammaPow))];
% contrast = [1 -1 0 0];
design = [gammaPow; ones(size(gammaPow))];

cfg=[];
% cfg.glm.contrast = contrast;
cfg.glm.statistic = 'beta';

for k=1:length(data_conc.label)
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

