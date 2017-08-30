function erf_osc_analysis_stat_glm_gamma(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time.mat', subj));
    tstat1{subj} = tmp{subj}.tstat1;
end
clear tmp

%% prepare for statistics

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
tstat1_GA = ft_appendtimelock(cfg, tstat1{allsubs});
clear tstat1

% create zero distribution
cfg = [];
cfg.parameter = 'trial';
cfg.operation = 'x1*0';
zero_distribution = ft_math(cfg, tstat1_GA);

%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
% neighbours           = ft_prepare_neighbours(cfg, tstat1_plCmb_GA); % define neighbouring channels
neighbours           = ft_prepare_neighbours(cfg, tstat1_GA); % define neighbouring channels

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

tstat2 = ft_timelockstatistics(cfg, tstat1_GA, zero_distribution);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_gamma_time_%s', erfoi);
save(fullfile([filename, '.mat']), 'tstat2', 'tstat1_GA');

ft_diary('off')



