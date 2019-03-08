function erf_osc_analysis_glm_gamma_virtualchan(subj, erfoi, filtdir)

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
ft_diary('on')
erf_osc_datainfo;

lcmvData = erf_osc_analysis_lcmv_orientation(subj, erfoi, filtdir, 'lp'); % optimal dipole orientation is chosen based on [0 0.5] after stimulus-reversal.
% data is timelocked to reversal or behavioral response, and low pass
% filtered. 
load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel.mat', subj, subj),'gammaPow'); % load gamma power

cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow=[-0.5 0];
lcmvData = ft_preprocessing(cfg, lcmvData);

if strcmp(erfoi, 'motor')
    cfg = [];
    cfg.latency = [-0.5 0];
    active = ft_selectdata(cfg, lcmvData);
else
    cfg = [];
    cfg.latency = [0 0.5];
    active = ft_selectdata(cfg, lcmvData);
    cfg.latency = [-0.5 0];
    baseline = ft_selectdata(cfg, lcmvData);
end

active.trial = cat(3,active.trial{:});
active.time = active.time{1};
if ~strcmp(erfoi, 'motor')
    baseline.trial = cat(3,baseline.trial{:});
    baseline.time = active.time;
end

% design = [gammaPow;((lcmvData.trialinfo(:,5)-lcmvData.trialinfo(:,4))/1200)'];
design=gammaPow;
cfg=[];
cfg.glm.statistic = 'beta';
cfg.glm.standardise = false;

dat = squeeze(active.trial);
dat = (dat - repmat(mean(dat,2),[1 length(lcmvData.trialinfo)]));
tmp = statfun_glm(cfg, dat, design);
betas_tmp = tmp.stat(:,1);

betas = rmfield(lcmvData, {'trial', 'time', 'trialinfo', 'sampleinfo'});
betas.time = active.time;
betas.avg = betas_tmp';
betas.dimord = 'chan_time';

if ~strcmp(erfoi, 'motor')
    dat_bl = squeeze(baseline.trial);
    dat_bl = (dat_bl - repmat(mean(dat_bl,2),[1 length(lcmvData.trialinfo)]));
    tmp_bl = statfun_glm(cfg, dat_bl, design);
    betas_bl_tmp = tmp_bl.stat(:,1);
    
    betas_bl = rmfield(betas, 'avg');
    betas_bl.avg = betas_bl_tmp';
else
    betas_bl=[];
end


filename = sprintf('/project/3011085.02/analysis/glm/sub-%03d/sub-%03d_glm_gamma_virtualchan_%s', subj, subj, erfoi);
save(fullfile([filename '.mat']), 'betas','betas_bl', '-v7.3');
ft_diary('off')

