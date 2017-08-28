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

%% regress out headmotion
load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/headmotion.mat', subj));
cfg=[];
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                              'HLC0021','HLC0022','HLC0023', ...
                              'HLC0031','HLC0032','HLC0033'};
headmotion = ft_selectdata(cfg, headmotion);

% calculate the mean coil position per trial
for trl = 1:nTrials
coil1(:,trl) = [mean(headmotion.trial{1,trl}(1,:)); mean(headmotion.trial{1,trl}(2,:)); mean(headmotion.trial{1,trl}(3,:))];
coil2(:,trl) = [mean(headmotion.trial{1,trl}(4,:)); mean(headmotion.trial{1,trl}(5,:)); mean(headmotion.trial{1,trl}(6,:))];
coil3(:,trl) = [mean(headmotion.trial{1,trl}(7,:)); mean(headmotion.trial{1,trl}(8,:)); mean(headmotion.trial{1,trl}(9,:))];
end
 
% calculate the headposition and orientation per trial (for function see bottom page) 
cc = circumcenter(coil1, coil2, coil3);

% demean to obtain translations and rotations from the average position and orientation
cc_dem = [cc - repmat(mean(cc,2),1,size(cc,2))]';


% add head movements to the regressorlist. also add the constant (at the end; column 7)
confound = [cc_dem ones(size(cc_dem,1),1)];

%% GLM on all trials

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
cfg=[];
cfg.latency = [-0.5 0];
baseline = ft_selectdata(cfg, data_conc);
cfg.latency = [0 0.5];
active = ft_selectdata(cfg, data_conc);


design = [gammaPow zeros(size(gammaPow)); zeros(size(gammaPow)) gammaPow; ones(size(gammaPow)) ones(size(gammaPow)); ...
    0.5*ones(size(gammaPow)) -0.5*ones(size(gammaPow)); cc_dem' cc_dem'];
contrast = [1 -1 0 0 0 0 0 0 0 0];

cfg=[];
cfg.glm.contrast = contrast;
cfg.glm.statistic = 'T';

for k=1:length(data_conc.label)
    dat = [squeeze(active.trial(:,k,:)); squeeze(baseline.trial(:,k,:))]';
    dat = (dat - repmat(mean(dat,2),[1 length(data.trialinfo)*2]))./(repmat(std(dat,[],2),[1 length(data.trialinfo)*2]));
    tmp = statfun_glm(cfg, dat, design);
    tstat1_tmp(k,:) = tmp.stat;
end

%% planar gradiant transformation of beta weights
% put beta weights in timelock structure
tstat1        = rmfield(data,{'trial', 'cfg'});
tstat1.avg    = tstat1_tmp;
tstat1.time   = active.time;
tstat1.dimord = 'chan_time';

% planar combination
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tstat1);
cfg.planarmethod    = 'sincos';
tstat1_planar            = ft_megplanar(cfg, tstat1);
cfg                 = [];
cfg.demean          = 'no';
tstat1_plCmb        = ft_combineplanar(cfg,tstat1_planar);


%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_gamma_time', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time', subj);
end
save(fullfile([filename '.mat']), 'tstat1', 'tstat1_plCmb', '-v7.3');
ft_diary('off')

