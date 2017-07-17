function erf_osc_analysis_stat_glm_erf(freqRange, zeropoint, erfoi)

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
if nargin<3
    erfoi = 'onset';
end
if isempty(erfoi)
    erfoi = 'onset';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_%s_erf_%s.mat', subj, freqRange, zeropoint, erfoi));
    avg_shuffles{subj} = tmp{subj}.shufflesAvgPlCmb;
    std_shuffles{subj} = tmp{subj}.shufflesStdPlCmb;
    betasPlCmb{subj} = tmp{subj}.betasPlCmb;
end

%% Baseline correct

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract'; %relative change
for subj=allsubs
    diffBetasPlCmb{subj} = ft_math(cfg, betasPlCmb{subj}, avg_shuffles{subj});
end

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
diffBetasPlCmbAvg = ft_appendfreq(cfg, diffBetasPlCmb{allsubs});
avg_betas_GA = ft_appendfreq(cfg, betasPlCmb{allsubs});
avg_shuffles_GA = ft_appendfreq(cfg, avg_shuffles{allsubs});



%% Do statistics

Nsub = length(allsubs);

cfg             = [];
cfg.method      = 'template'; 
cfg.feedback    = 'no';
neighbours      = ft_prepare_neighbours(cfg, avg_betas_GA); % define neighbouring channels

cfg                  = [];
cfg.channel          = 'MEG';
cfg.neighbours       = neighbours;
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
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

stat = ft_freqstatistics(cfg, avg_betas_GA, avg_shuffles_GA);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_tf_%s_%s_erf_%s.mat', freqRange, zeropoint, erfoi);
save(filename, 'stat', 'diffBetasPlCmbAvg');

ft_diary('off')



