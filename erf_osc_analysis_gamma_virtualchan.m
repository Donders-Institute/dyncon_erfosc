 function erf_osc_analysis_gamma_virtualchan(subj, isPilot, sourcemodel)
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
    sourcemodel = '3d'; % note that this variable is replaced by the actual sourcemodel
end
if isempty(sourcemodel)
    sourcemodel = '3d';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(fullfile([pilotsubjects(subj).segmentedmri, '.mat']));
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak', subj), 'peakFreq_gamma');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(fullfile([subjects(subj).mridir, 'preproc/headmodel.mat']));
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj), 'peakFreq_gamma');
    if strcmp(sourcemodel, '2d')
        load(fullfile([subjects(subj).mridir, 'preproc/sourcemodel2d.mat']));
    else
        load(fullfile([subjects(subj).mridir, 'preproc/sourcemodel3d.mat']));
    end
end

data = data.dataClean;
fs = data.fsample;

cfg         = [];
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM;
data       = ft_selectdata(cfg, data);

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4)); % trialinfo is specified in 1200 Hz. If data is resampled, it has to be taken care of for ft_redefinetrial.
dataShift  = ft_redefinetrial(cfg, data);


% select data: 0.25 second preceding grating start and 0.25 second preceding
% grating shift. Contrast these in terms of gamma frequency at gamma peak
cfg         = [];
cfg.latency = [-0.25+1/fs 0];
dataPre     = ft_selectdata(cfg, data);
% take second preceding shift (NOTE: is it confounding if this includes
dataPost    = ft_selectdata(cfg, dataShift);
dataAll     = ft_appenddata([], dataPost, dataPre);

%% Frequency analysis
% calculate power in complete data (for calculation of common filter)
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 8;
cfg.foilim    = [peakFreq_gamma peakFreq_gamma];
freqAll       = ft_freqanalysis(cfg, dataAll);

% calculate power pre and post stimulus
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'powandcsd';
cfg.tapsmofrq  = 8;
cfg.foilim     = [peakFreq_gamma peakFreq_gamma];
freqPre        = ft_freqanalysis(cfg, dataPre);
cfg.keeptrials = 'yes';
% cfg.pad        = 6;
freqPost       = ft_freqanalysis(cfg, dataPost);
% cfg.tapsmofrq  = 1;
% cfg.foilim     = [peakFreq_alpha peakFreq_alpha];
% freqPost_alpha = ft_freqanalysis(cfg, dataPost);

%% Source analysis

% computes the forward model for many dipole locations on a regular 2D or
% 3D grid and stores it for efficient inverse modelling
cfg                 = [];
cfg.grad            = data.grad;
cfg.headmodel       = headmodel;
cfg.reducerank      = 2;
cfg.grid            = sourcemodel; % use the grid specified in sourcemodel
cfg.channel         = 'MEG';
leadfield           = ft_prepare_leadfield(cfg, freqAll);

% calculate common filter
cfg                   = [];
cfg.grad              = freqAll.grad;
cfg.method            = 'dics';
cfg.frequency         = peakFreq_gamma;
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';
sourceAll             = ft_sourceanalysis(cfg, freqAll);

cfg.grid.filter       = sourceAll.avg.filter;
sourcePre             = ft_sourceanalysis(cfg, freqPre);
sourcePost            = ft_sourceanalysis(cfg, freqPost);
sourceDiff            = sourcePost;
sourceDiff.avg.pow    = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;

%% Create virtual channel of maximal gamma power
% the following position contains the max gamma power difference
% [maxval, maxpowindx] = max(sourceDiff.avg.pow(sourceDiff.inside));
[maxval, maxpowindx] = max(sourceDiff.avg.pow);
sourceDiff.pos(maxpowindx, :)
% we will create a virtual channel based on this location. In order to do
% this, we have to do timelockanalysis and use an LCMV beamformer, because
% this will pass the activity at the location of interest with unit gain,
% while optimally reducing activity from other locations. This filter then
% can be applied to the original MEG data

% select the grid information only for the position with maximum gamma
% power
virtualgrid = cfg.grid;
virtualgrid.pos = virtualgrid.pos(maxpowindx,:);
virtualgrid.inside = virtualgrid.inside(maxpowindx);
virtualgrid.leadfield = {[virtualgrid.leadfield{maxpowindx}]};
virtualgrid.filter = {[virtualgrid.filter{maxpowindx}]};

cfg.grid = virtualgrid;
cfg.rawtrial='yes';
gammaChan = ft_sourceanalysis(cfg, freqPost);

% alphaChan = ft_sourceanalysis(cfg, freqPost_alpha);

cfg                   = [];
cfg.covariance        = 'yes';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, dataAll);

cfg                 = [];
cfg.method          = 'lcmv';
cfg.headmodel       = headmodel;
cfg.grid.pos        = sourcemodel.pos(maxpowindx,:);
cfg.grid.inside     = sourcemodel.inside(maxpowindx);
cfg.grid.unit       = sourcemodel.unit;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectmom = 'yes';
cfg.lcmv.fixedori   = 'yes';
source_idx          = ft_sourceanalysis(cfg, tlock);

gammaFilter = source_idx.avg.filter;

lcmvData              = data;
lcmvData.label        = {'gam_pow'};
lcmvData.trial        = [];
for i=1:length(dataShift.trial)
    lcmvData.trial{i} = gammaFilter{1} * data.trial{i};
end

%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel', subj);
end
save(fullfile([filename '.mat']), 'lcmvData', 'gammaFilter', 'gammaChan');
ft_diary('off')


end

