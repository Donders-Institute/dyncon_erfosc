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
    avg_betas{subj} = tmp{subj}.avg_shuffles_struct;
    std_betas{subj} = tmp{subj}.std_shuffles_struct;
end

%% Baseline correct
for subj=allsubs
    if strcmp(freqRange, 'high')
        betasPlCmb{subj} = ft_freqdescriptives([], tmp{subj}.bhPlanarCmb);
    else
        betasPlCmb{subj} = ft_freqdescriptives([], tmp{subj}.blPlanarCmb);
    end
end

% "z-score" with shuffled beta weights.
cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)./x3'; %relative change
for subj=allsubs
    diffBetasPlCmb{subj} = ft_math(cfg, betasPlCmb{subj}, avg_betas{subj}, std_betas{subj});
end

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
diffBetasPlCmbAvg = ft_appendfreq(cfg, diffBetasPlCmb{allsubs});

zerodistribution = diffBetasPlCmbAvg;
zerodistribution.powspctrm = zerodistribution.powspctrm*0;


%% Do statistics

Nsub = length(allsubs);

cfg             = [];
cfg.method      = 'template'; 
cfg.feedback    = 'no';
neighbours      = ft_prepare_neighbours(cfg, diffBetasPlCmbAvg); % define neighbouring channels

cfg                  = [];
cfg.channel          = 'MEG';
cfg.neighbours       = neighbours;
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
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

stat = ft_freqstatistics(cfg, diffBetasPlCmbAvg, zerodistribution);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_tf_%s_%s_erf_%s.mat', freqRange, zeropoint, erfoi);
save(filename, 'stat', 'diffBetasPlCmbAvg');

ft_diary('off')



