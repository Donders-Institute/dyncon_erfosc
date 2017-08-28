function erf_osc_analysis_glm_gamma_time(subj, isPilot, erfoi, doDSS)
% do a linear regression of pre-change gamma power over time.
if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(erfoi)
    erfoi = 'reversal';
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
    end
end

for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
end
gammaPow = (gammaPow-mean(gammaPow))/std(gammaPow);
nTrials = length(data.trial);

[~, idxMax] = sort(gammaPow, 2, 'descend');

%% GLM on all trials
design = [ones(size(gammaPow)); gammaPow];
cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 30;
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

for k=1:length(data_conc.label)
    Y = squeeze(data_conc.trial(:,k,:));
    betas(:,:,k) = design'\Y;
end

%% planar gradiant transformation of beta weights
% put beta weights in timelock structure
tl=[];
tl.avg    = squeeze(betas(2,:,:))';
tl.time   = data.time{1};
tl.dimord = 'chan_time';
tl.label  = data.label;
tl.grad   = data.grad;

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
t1        = nearest(tlPlanarCmb.time, -0.25);
t2        = nearest(tlPlanarCmb.time, 0);
mu        = rmfield(tlPlanarCmb, 'avg');
mu.avg    = mean(tlPlanarCmb.avg(:,t1:t2), 2);
mu.avg    = repmat(mu.avg, [1, length(tlPlanarCmb.time)]);
sigma     = rmfield(tlPlanarCmb, 'avg');
sigma.avg = std(tlPlanarCmb.avg(:,t1:t2),[],2);
sigma.avg = repmat(sigma.avg, [1, length(tlPlanarCmb.time)]);

cfg=[];
cfg.parameter = 'avg';
cfg.operation = '(x1-x2)./x3';
tlPlanarCmbZ = ft_math(cfg, tlPlanarCmb, mu, sigma);

%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_gamma_time', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time', subj);
end
save(fullfile([filename '.mat']), 'betas', 'tlPlanarCmb','tlPlanarCmbZ','tl', '-v7.3');
ft_diary('off')

