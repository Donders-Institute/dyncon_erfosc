cfg                 = [];
cfg.method          = 'broadband';
cfg.output          = 'all';
cfg.randomseed      = 'yes';
cfg.numtrl          = 500;
cfg.trllen          = 1;
cfg.fsample         = 1200;

cfg.n1.ampl         = 1;
cfg.n1.bpfreq       = [5 100];


cfg.noise.ampl      = 0.25;

bl = ft_freqsimulation(cfg);

cfg.n2.ampl         = 0.1;
cfg.n2.bpfreq       = [40 60];

act = ft_freqsimulation(cfg);

%% freqanalysis
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foilim      = [2 100];
cfg.keeptrials  = 'no'; % average baseline over trials
cfg.channel     = 'mix';
pow_bl          = ft_freqanalysis(cfg, bl);
pow_act         = ft_freqanalysis(cfg, act);