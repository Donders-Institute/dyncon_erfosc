function erf_osc_analysis_stat_tfr(freqRange, zeropoint)

if nargin<1
    freqRange = 'high';
end
if isempty(freqRange)
    freqRange = 'high';
end
if nargin<2
    zeropoint = 'onset';
end
if isempty(zeropoint)
    zeropoint = 'onset';
end
ft_diary('on')

erf_osc_datainfo;

% load data
for subj=allsubs
    if strcmp(freqRange, 'high')
        tfr{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s.mat', subj, zeropoint), 'tfaHigh');
        tfr{subj} = tfr{subj}.tfaHigh;
    else
        tfr{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s.mat', subj, zeropoint), 'tfaLow');
        tfr{subj} = tfr{subj}.tfaLow;
    end
    bl{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_onset.mat', subj), 'baselineH');
    bl{subj} = bl{subj}.baselineH;
end

for subj=allsubs
    tfr{subj} = ft_freqdescriptives([], tfr{subj});
    bl{subj}.time = tfr{subj}.time;
    bl{subj}.powspctrm = repmat(bl{subj}.powspctrm, [1,1,length(bl{subj}.time)]);
end

cfg=[];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract';
for subj=allsubs
    diff{subj} = ft_math(cfg, tfr{subj}, bl{subj});
end

diffAvg = ft_freqgrandaverage([], diff{allsubs});
cfg=[];
cfg.keepindividual = 'yes';
tfrAvg = ft_freqgrandaverage(cfg, tfr{allsubs});
blAvg = ft_freqgrandaverage(cfg, tfr{allsubs});


%% statistics

Nsub = length(allsubs);

cfg             = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.feedback    = 'no';
neighbours      = ft_prepare_neighbours(cfg, tfrAvg); % define neighbouring channels

cfg                  = [];
cfg.channel          = 'MEG';
cfg.neighbours       = neighbours;
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_actvsblT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.numrandomization = 10000;

% cfg.design
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg, tfrAvg, blAvg);


% save
filename = '/project/3011085.02/results/stat_tfr.mat';
save(filename, 'stat');

ft_diary('off')

