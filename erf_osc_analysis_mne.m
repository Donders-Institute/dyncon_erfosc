function erf_osc_analysis_mne(subj, isPilot, sourcemodel)
% This function estimates gamma power at the estimated gamma peak frequency

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
    sourcemodel = '2d'; % note that this variable is replaced by the actual sourcemodel
end
if isempty(sourcemodel)
    sourcemodel = '2d';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
    load(fullfile([pilotsubjects(subj).segmentedmri, '.mat']));
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
    load(fullfile([subjects(subj).mridir, 'preproc/headmodel.mat']));
    if strcmp(sourcemodel, '2d')
        load(fullfile([subjects(subj).mridir, 'preproc/sourcemodel2d.mat']));
    else
        load(fullfile([subjects(subj).mridir, 'preproc/sourcemodel3d.mat']));
    end
    sourcemodel = ft_convert_units(sourcemodel, 'mm');
end

data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
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

% preprocess data for lmcv filter
% lowpass filter stimulus reversal erf
cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.lpfilttype = 'firws';
data = ft_preprocessing(cfg, data);

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift  = ft_redefinetrial(cfg, data);
cfg.offset = -(data.trialinfo(:,6)-data.trialinfo(:,4)); % trialinfo is specified in 1200 Hz. If data is resampled, it has to be taken care of for ft_redefinetrial.
dataResp  = ft_redefinetrial(cfg, data);

% select data: 0.75 second preceding grating start same for grating shift.
% Contrast these in terms of gamma frequency at gamma peak
cfg         = [];
cfg.latency = [-0.2+1/fs 0];
dataPreResp = ft_selectdata(cfg, dataResp);
dataPreStim = ft_selectdata(cfg, data); % use this for finding source location
dataAll     = ft_appenddata([], dataPreResp, dataPreStim);

cfg=[];
cfg.vartrllength = 2;
cfg.keeptrials = 'yes';
tlckPreResp = ft_timelockanalysis(cfg, dataPreResp)
tlckAll     = ft_timelockanalysis(cfg, dataAll);
cfg.covariance = 'yes';
tlckPreStim = ft_timelockanalysis(cfg, dataPreStim);

tlckPreResp.cov = tlckPreStim.cov;
%% Compute covariance matrix (over toi and baseline together)
cfg        = [];
cfg.method = 'mne';
cfg.grid   = sourcemodel;
cfg.headmodel    = headmodel;
cfg.mne.prewhiten = 'yes';
cfg.mne.lambda    = 3;
cfg.mne.scalesourcecov = 'yes';
sourcePreResp  = ft_sourceanalysis(cfg,tlckPreResp);
sourcePreStim = ft_sourceanalysis(cfg, tlckPreStim);


% compute spatial filter based on entire data
cfg = [];
cfg.covariance='yes';
cfg.covariancewindow = [-0.2+1/fs 0];
cfg.vartrllength = 2;
avg = ft_timelockanalysis(cfg,dataAll);

cfg = [];
cfg.covariance='yes';
avgresp = ft_timelockanalysis(cfg,dataPreResp);
avgbl = ft_timelockanalysis(cfg,dataPreStim);

%% source analysis
% compute spatial filter
cfg                 = [];
cfg.method          = 'lcmv';
cfg.grid            = sourcemodel;
cfg.headmodel       = headmodel;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda     = '5%';
cfg.channel         = {'MEG'};
cfg.senstype        = 'MEG';
sourceavg           = ft_sourceanalysis(cfg, avg);

% second and third call to ft_sourceanalysis now applying the precomputed filters to pre and %post intervals
cfg             = [];
cfg.method      = 'lcmv';
cfg.grid        = sourcemodel;
cfg.grid.filter = sourceavg.avg.filter;
cfg.vol         = headmodel;
sourceResp      = ft_sourceanalysis(cfg, avgresp);
sourceBl        = ft_sourceanalysis(cfg, avgbl);

M1=sourceResp;
M1.avg.pow=(sourceResp.avg.pow-sourceBl.avg.pow)./sourceBl.avg.pow;
%% save
if isPilot
    filename = sprintf('/project/3011085.02/analysis/freq/pilot-%03d/sub-%03d_gamma_virtual_channel', subj, subj);
else
    filename = sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel', subj, subj);
end
save(fullfile([filename '.mat']), 'lcmvData', 'gammaFilter', 'gammaPow', 'sourceDiff', 'maxpowindx');
ft_diary('off')


end

