function erf_osc_analysis_stat_glm_gamma(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/analysis/GLM/sub-%03d/sub-%03d_glm_gamma_time_%s.mat', subj, subj, erfoi));
    betas{subj} = rmfield(tmp{subj}.betas, 'elec');
end
clear tmp

%% prepare for statistics

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
betas_GA = ft_appendtimelock(cfg, betas{allsubs});

% contrast with zero
ref = betas_GA;
ref.trial = ref.trial*0; % since there is no biased estimate of the regression coefficients, we can compare with zero instead of with baseline.

%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
neighbours           = ft_prepare_neighbours(cfg, betas_GA); % define neighbouring channels

cfg                  = [];
cfg.channel          = 'MEG';
cfg.neighbours       = neighbours;
cfg.parameter        = 'trial';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.numrandomization = 10000;


cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, betas_GA, ref);
stat.cfg = rmfield(stat.cfg, 'previous');


% save
filename = sprintf('/project/3011085.02/analysis/stat_glm_gamma_time_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat');

ft_diary('off')



