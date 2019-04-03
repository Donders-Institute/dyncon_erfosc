function erfosc_analysis_gamma_virtualchan(subj, isPilot, sourcemodel)
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
erfosc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
    load(fullfile([pilotsubjects(subj).segmentedmri, '.mat']));
    load(sprintf('/project/3011085.02/analysis/freq/pilot-%03d/sub-%03d_gamma_peak', subj, subj), 'peakFreq_gamma');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
    load(fullfile([subjects(subj).mridir, 'preproc/headmodel.mat']));
    load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_pow.mat', subj, subj), 'peakFreq_gamma');
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

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4)); % trialinfo is specified in 1200 Hz. If data is resampled, it has to be taken care of for ft_redefinetrial.
dataShift  = ft_redefinetrial(cfg, data);


% select data: 0.75 second preceding grating start same for grating shift.
% Contrast these in terms of gamma frequency at gamma peak
cfg         = [];
cfg.latency = [-0.75+1/fs 0];
dataPreStim = ft_selectdata(cfg, data);
dataPreRev  = ft_selectdata(cfg, dataShift); % use this for finding source location
cfg.latency = [-0.25+1/fs 0];
dataPreRev_short = ft_selectdata(cfg, dataShift);
dataAll     = ft_appenddata([], dataPreRev, dataPreStim);


% preprocess data for lmcv filter
% lowpass filter stimulus reversal erf
cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.lpfilttype = 'firws';
dataPostRev = ft_preprocessing(cfg, dataShift);

cfg=[];
cfg.latency = [-0.2+1/fs 0.5];
dataPostRev = ft_selectdata(cfg, dataPostRev);
%% Frequency analysis
% calculate power in complete data (for calculation of common filter)
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 8;
cfg.pad       = 1;
cfg.foilim    = [peakFreq_gamma peakFreq_gamma];
freqAll       = ft_freqanalysis(cfg, dataAll);

% calculate power pre and post stimulus
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'fourier';%'powandcsd';
cfg.tapsmofrq    = 8;
cfg.pad          = 1;
cfg.foilim       = [peakFreq_gamma peakFreq_gamma];
freqPreStim      = ft_freqanalysis(cfg, dataPreStim);
cfg.keeptrials   = 'yes';
freqPreRev       = ft_freqanalysis(cfg, dataPreRev);
freqPreRev_short = ft_freqanalysis(cfg, dataPreRev_short);



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
sourcePre             = ft_sourceanalysis(cfg, freqPreStim);
sourcePost            = ft_sourceanalysis(cfg, freqPreRev);
sourceDiff            = sourcePost;
sourceDiff.avg.pow    = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;

%% Create virtual channel of maximal gamma power
% the following position contains the max gamma power difference
% [maxval, maxpowindx] = max(sourceDiff.avg.pow(sourceDiff.inside));
[maxval, maxpowindx] = max(sourceDiff.avg.pow);
if subj==15
    [val, idx] = sort(sourceDiff.avg.pow, 'descend');
    idx(isnan(val))=[];
    val(isnan(val))=[];
    maxpowindx = idx(39);
end
sourceDiff.pos(maxpowindx, :)
% we will create a virtual channel based on this location. In order to do
% this, we have to do timelockanalysis and use an LCMV beamformer, because
% this will pass the activity at the location of interest with unit gain,
% while optimally reducing activity from other locations. This filter then
% can be applied to the original MEG data
%
% select the grid information only for the position with maximum gamma
% power
virtualgrid = rmfield(cfg.grid, 'filter');
virtualgrid.pos = virtualgrid.pos(maxpowindx,:);
virtualgrid.inside = virtualgrid.inside(maxpowindx);
virtualgrid.leadfield = {[virtualgrid.leadfield{maxpowindx}]};


cfg.method = 'pcc';
cfg.pcc = cfg.dics;
cfg = rmfield(cfg, 'dics');
cfg.grid = virtualgrid;
%cfg.rawtrial='yes';
gammaChan = ft_sourceanalysis(cfg, freqPreRev_short); % should the filter just be determined on (250ms or 750ms) gamma data?
% here we have the single taper fourier coefficients in the mom-field
mom = gammaChan.avg.mom{1};
[u,s,v]=svd(real(mom*mom'));
newmom=u(:,1)'*mom;
newpow=abs(newmom).^2;
newpow_trl=[];
nTapers = size(mom,2)/length(dataPreRev_short.trial);
for iTrial = 1:length(dataPreRev_short.trial)
    newpow_trl(iTrial) = sum(newpow(iTrial*nTapers-(nTapers-1):iTrial*nTapers))/nTapers;
end
gammaPow = log(newpow_trl);
gammaPow = (gammaPow-mean(gammaPow))/std(gammaPow);


%%%%%% Previously %%%%%%%
%{
virtualgrid = cfg.grid;
virtualgrid.pos = virtualgrid.pos(maxpowindx,:);
virtualgrid.inside = virtualgrid.inside(maxpowindx);
virtualgrid.leadfield = {[virtualgrid.leadfield{maxpowindx}]};
virtualgrid.filter = {[virtualgrid.filter{maxpowindx}]};

cfg.grid = virtualgrid;
cfg.rawtrial='yes';
gammaChan = ft_sourceanalysis(cfg, freqPreRev_short);
for k=1:length(gammaChan.trial)
    gammaPow(k) = gammaChan.trial(k).pow;
end
gammaPow = log(gammaPow);
gammaPow = gammaPow-mean(gammaPow);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%


cfg                   = [];
cfg.covariance        = 'yes';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, dataPostRev);

cfg                 = [];
cfg.method          = 'lcmv';
cfg.headmodel       = headmodel;
cfg.grid.pos        = sourcemodel.pos(maxpowindx,:);
cfg.grid.inside     = sourcemodel.inside(maxpowindx);
cfg.grid.unit       = sourcemodel.unit;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda     = '5%';
source_idx          = ft_sourceanalysis(cfg, tlock);

gammaFilter = source_idx.avg.filter;

lcmvData              = data;
lcmvData.label        = {'gam_pow_1', 'gam_pow_2', 'gam_pow_3'};
lcmvData.trial        = [];
for i=1:length(dataShift.trial)
    lcmvData.trial{i} = gammaFilter{1} * data.trial{i};
end


sourceDiff = rmfield(sourceDiff,'cfg');

%% save
if isPilot
    filename = sprintf('/project/3011085.02/analysis/freq/pilot-%03d/sub-%03d_gamma_virtual_channel', subj, subj);
else
    filename = sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel', subj, subj);
end
save(fullfile([filename '.mat']), 'lcmvData', 'gammaFilter', 'gammaPow', 'sourceDiff', 'maxpowindx');

ft_diary('off')


end

