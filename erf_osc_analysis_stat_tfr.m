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
cfg.operation = '(x1-x2)./x2';
for subj=allsubs
    diff{subj} = ft_math(cfg, tfr{subj}, bl{subj});
end

cfg = [];
cfg.appenddim = 'rpt';
tfrAvg        = ft_appendfreq(cfg, tfr{allsubs});
blAvg         = ft_appendfreq(cfg, bl{allsubs});

cfg=[];
diffAvg       = ft_freqgrandaverage(cfg, diff{allsubs});


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
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.01;
cfg.correcttail      = 'prob';
cfg.minnbchan        = 2;
cfg.numrandomization = 10000;

% cfg.design
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg, tfrAvg, blAvg);


% save
filename = sprintf('/project/3011085.02/results/stat_tfr_%s.mat', zeropoint);
save(filename, 'stat', 'diffAvg', '-v7.3');

ft_diary('off')

