function erf_osc_analysis_stat_glm_erf(freqRange, zeropoint, erfoi)

if nargin<1 || isempty(freqRange)
    freqRange = 'high';
end
if nargin<2 || isempty(zeropoint)
    zeropoint = 'reversal';
end
if nargin<3 || isempty(erfoi)
    erfoi = 'reversal';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_%s_erf_%s.mat', subj, freqRange, zeropoint, erfoi));
    avg_shuffles{subj} = tmp{subj}.shufflesAvg;
    std_shuffles{subj} = tmp{subj}.shufflesStd;
    betas{subj} = tmp{subj}.betas;
    betas_baseline{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_onset_erf_%s.mat',subj, freqRange, erfoi),'betas_baseline');
    betas_baseline{subj} = betas_baseline{subj}.betas_baseline;
    betas_baseline_no_avg = betas_baseline;
end
clear tmp
%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVE VS. BASELINE %
%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline correct
% average baseline over time and repeat over time

for subj = allsubs
    cfg=[];
    cfg.avgovertime = 'yes';
    betas_baseline{subj} = ft_selectdata(cfg, betas_baseline{subj});
    betas_baseline{subj}.powspctrm = repmat(betas_baseline{subj}.powspctrm, [1,1, length(betas{1}.time)]);
    betas_baseline{subj}.time = betas{subj}.time;
    cfg=[];
    cfg.latency=[betas{1}.time(end-15) betas{1}.time(end)];
    betas_short{subj} = ft_selectdata(cfg, betas{subj});
    betas_baseline_no_avg{subj}.time = betas_short{subj}.time;
end

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract'; %relative change
for subj=allsubs
    diffBetas{subj} = ft_math(cfg, betas{subj}, betas_baseline{subj});
    diffBetas_no_avg{subj} = ft_math(cfg, betas_short{subj}, betas_baseline_no_avg{subj});
end

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
diffBetasAvg = ft_appendfreq(cfg, diffBetas{allsubs});
diffBetasAvg2 = ft_appendfreq(cfg, diffBetas_no_avg{allsubs});
betas_GA = ft_appendfreq(cfg, betas{allsubs});
betas_short_GA = ft_appendfreq(cfg, betas_short{allsubs});
betas_baseline_GA = ft_appendfreq(cfg, betas_baseline{allsubs});
betas_baseline_no_avg_GA = ft_appendfreq(cfg, betas_baseline_no_avg{allsubs});

%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVE VS. SHUFFLES %
%%%%%%%%%%%%%%%%%%%%%%%
%{
%% Baseline correct

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract'; %relative change
for subj=allsubs
    diffBetas{subj} = ft_math(cfg, betas{subj}, avg_shuffles{subj});
end

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
diffBetasAvg = ft_appendfreq(cfg, diffBetas{allsubs});
avg_betas_GA = ft_appendfreq(cfg, betas{allsubs});
avg_shuffles_GA = ft_appendfreq(cfg, avg_shuffles{allsubs});

%}

%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
neighbours           = ft_prepare_neighbours(cfg, betas_GA); % define neighbouring channels

cfg                  = [];
cfg.channel          = 'MEG';
cfg.neighbours       = neighbours;
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.minnbchan        = 2;
cfg.correcttail      = 'prob';
cfg.numrandomization = 10000;


cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg, betas_GA, betas_baseline_GA);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_tf_%s_%s_erf_%s_2.mat', freqRange, zeropoint, erfoi);
save(filename, 'stat', 'diffBetasAvg', 'betas_GA', 'betas_baseline_GA', 'betas_baseline_no_avg_GA', 'betas_short_GA');

ft_diary('off')



