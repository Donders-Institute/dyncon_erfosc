function erf_osc_preprocessing_ica(subj, isPilot)
% This function decomposes the data with ICA and should be called before
% erf_osc_preprocessing_artifact. 50 components are computed with fastica
% procedure and saved on harddisk.

if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

%% Load data and define trials
erf_osc_datainfo; % load subject specific info.

cfg=[];
if isPilot
    cfg.dataset = pilotsubjects(subj).dataset;
    cfg.logfile = load(pilotsubjects(subj).logfile);% load log file
else
    cfg.dataset = subjects(subj).dataset;
    cfg.logfile = load(subjects(subj).logfile);
end
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erf_osc_mytrialfun';
cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);

% preprocess data
cfg.continuous = 'yes';
cfg.demean = 'yes';
cfg.padding = 8; % pad the data so the filter is the same for all trials (filter depends on trial length)
cfg.dftfilter = 'yes';
cfg.dftfreq   = [50+(-1:1)./cfg.padding 100 150 200]; % remove line noise
cfg.usefftfilt = 'yes'; % frequency multiplication of filter instead of convolution (slow)
cfg.hpfilter = 'yes';
cfg.hpfreq = 1; % remove slow drifts
cfg.hpfilttype = 'firws';
data = ft_preprocessing(cfg);

cfg=[];
cfg.resamplefs = 200;
cfg.detrend = 'no';
dataResample = ft_resampledata(cfg, data);

%% perform the independent component analysis (i.e., decompose the data)
cfg=[];
cfg.method          = 'fastica';
cfg.channel         = 'MEG';
cfg.fastica.numOfIC = 50;
comp = ft_componentanalysis(cfg, dataResample);

% save
if isPilot
    save(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_icaComp.mat', subj,subj), 'comp');
else
    save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_icaComp.mat', subj, subj), 'comp');
end

end

