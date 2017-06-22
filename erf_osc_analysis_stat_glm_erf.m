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
    if strcmp(zeropoint, 'reversal')
        bl{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_onset_erf_%s.mat', subj, freqRange, erfoi));
        bl{subj} = bl{subj}.bhPlanarCmb;
    end
end

%% Baseline correct
for subj=allsubs
    betasPlCmb{subj} = ft_freqdescriptives([], tmp{subj}.bhPlanarCmb);
end

% get active period
cfg             = [];
if strcmp(zeropoint, 'onset')
    cfg.latency = [0 1.5];
else
    cfg.latency = [-1 0.5];
end
for subj=allsubs
    betasPlCmb_act{subj} = ft_selectdata(cfg, betasPlCmb{subj});
end

% get baseline
% for data timelocked to stimulus reversal, take baseline in data locked to stimulus onset 
cfg             = [];
cfg.avgovertime = 'yes';
cfg.latency     = [-1 -0.2];

for subj=allsubs
    if strcmp(zeropoint, 'onset');
        betasPlCmb_bl{subj} = ft_selectdata(cfg, betasPlCmb{subj});
    else
        betasPlCmb_bl{subj} = ft_selectdata(cfg, bl{subj});
    end
    betasPlCmb_bl{subj}.time = betasPlCmb_act{subj}.time;
    betasPlCmb_bl{subj}.powspctrm = repmat(betasPlCmb_bl{subj}.powspctrm, [1, 1, length(betasPlCmb_bl{subj}.time)]);
end

% baseline correct
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)./x2'; %relative change
for subj=allsubs
    diffBetasPlCmb{subj} = ft_math(cfg, betasPlCmb_act{subj}, betasPlCmb_bl{subj});
end

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
betasPlCmbAvg_bl  = ft_appendfreq(cfg, betasPlCmb_bl{allsubs});
betasPlCmbAvg_act = ft_appendfreq(cfg, betasPlCmb_act{allsubs});
diffBetasPlCmbAvg = ft_appendfreq(cfg, diffBetasPlCmb{allsubs});


%% Do statistics

Nsub = length(allsubs);

cfg             = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.feedback    = 'no';
neighbours      = ft_prepare_neighbours(cfg, betasPlCmbAvg_act); % define neighbouring channels

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

stat = ft_freqstatistics(cfg, betasPlCmbAvg_act, betasPlCmbAvg_bl);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_tf_%s_%s_erf_%s.mat', freqRange, zeropoint, erfoi);
save(filename, 'stat', 'diffBetasPlCmbAvg');

ft_diary('off')



