function powBasic(subj)
subj=3;
%% powBasic
isPilot=true;

erf_osc_datainfo;
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


if subj==103
    cfg.trialdef.poststim = 2;
    cfg.trialdef.prestim = 1;
    cfg.trialfun = 'erf_osc_mytrialfun_gamma_stan';
else
    cfg.trialfun = 'erf_osc_mytrialfun2';
    cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
    cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
end
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
cfg.demean = 'yes';
data = ft_preprocessing(cfg);
%%
fs = data.fsample;
cfg=[];
cfg.latency = [-0.5+1/fs 0];
dataBaseline = ft_selectdata(cfg, data);
cfg.latency = [0.3+1/fs 0.8]; % take active time window after first erfs
dataActive  = ft_selectdata(cfg, data);

cfg=[];
cfg.foilim = [2 100];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.tapsmofrq = 8;
cfg.taper = 'dpss';
cfg.channel = 'MEG';
cfg.keeptrials = 'yes';
powActive = ft_freqanalysis(cfg, dataActive);
powBaseline = ft_freqanalysis(cfg, dataBaseline);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 6;
cfg.foi          = 2:2:30;% analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -1:0.05:3.75;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials   = 'yes';
tfaLow = ft_freqanalysis(cfg, data);


cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 8; % 8 Hz freq smoothing on both sides
cfg.foi          = 28:4:100;
cfg.pad          = 6;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*(1/4);
cfg.toi          = -1.5:0.05:3.75;
cfg.keeptrials   = 'yes';
tfaHigh = ft_freqanalysis(cfg,data);

if subj==103
    save('pilot3_tfa_stan', 'powActive', 'powBaseline', 'tfaLow', 'tfaHigh', '-v7.3');
else
    save('pilot3_tfa', 'powActive', 'powBaseline', 'tfaLow', 'tfaHigh', '-v7.3');
end
end
