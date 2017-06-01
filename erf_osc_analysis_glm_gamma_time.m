function erf_osc_analysis_glm_gamma_time(subj, isPilot)
% do a linear regression of pre-change gamma power over time.
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

% close all

% Initiate Diary
ft_diary('on')


%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/dss.mat', subj), 'data_dss');
    %     load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
end
fs=data_dss.fsample;
for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
%     gammaPow2(i) = (gammaChan.trial(i).pow);
end
gammaPow = gammaPow-mean(gammaPow);
nTrials = length(data_dss.trial);

[~, idxMax] = sort(gammaPow, 2, 'descend');

%% GLM on all trials
design = [ones(size(gammaPow)); gammaPow];
data=data_dss;
cfg=[];
cfg.lpifilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 40;
data = ft_preprocessing(cfg, data);
data.trial = cat(3,data.trial{:});
data.trial = permute(data.trial, [3,1,2]);
data.time = data.time{1};

for k=1:length(data.label)
    Y = squeeze(data.trial(:,k,:));
    betas(:,:,k) = design'\Y;
end

%% planar gradiant transformation of beta weights
% put beta weights in timelock structure
tl=[];
tl.avg    = squeeze(betas(2,:,:))';
tl.time   = data_dss.time{1};
tl.dimord = 'chan_time';
tl.label  = data_dss.label;
tl.grad   = data_dss.grad;

% planar combination
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
cfg.planarmethod    = 'sincos';
tlPlanar            = ft_megplanar(cfg, tl);
cfg                 = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.25 0];
tlPlanarCmb         = ft_combineplanar(cfg,tlPlanar);

%% Normalize beta weights
% zscore manually based on baseline window
t1        = nearest(tl.time, -0.25);
t2        = nearest(tl.time, 0);
mu        = rmfield(tl, 'avg');
mu.avg    = mean(tl.avg(:,t1:t2), 2);
mu.avg    = repmat(mu.avg, [1, length(tl.time)]);
sigma     = rmfield(tl, 'avg');
sigma.avg = std(tl.avg(:,t1:t2),[],2);
sigma.avg = repmat(sigma.avg, [1, length(tl.time)]);

cfg=[];
cfg.parameter = 'avg';
cfg.operation = '(x1-x2)./x3';
tlPlanarCmbZ = ft_math(cfg, tl, mu, sigma);

%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_gamma_time', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time', subj);
end
save(fullfile([filename '.mat']), 'betas', 'tlPlanarCmbZ', '-v7.3');
ft_diary('off')

