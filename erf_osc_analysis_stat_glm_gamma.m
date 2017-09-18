function erf_osc_analysis_stat_glm_gamma(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time_%s.mat', subj, erfoi));
    betas_plcmb{subj} = tmp{subj}.betas_plcmb;
    if strcmp(erfoi, 'motor')
        betas_bl_plcmb{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time_%s.mat', subj, 'reversal'), 'betas_bl_plcmb');
        betas_bl_plcmb{subj} = betas_bl_plcmb{subj}.betas_bl_plcmb;
    else
        betas_bl_plcmb{subj} = tmp{subj}.betas_bl_plcmb;
    end
end
clear tmp

%% prepare for statistics

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
betas_plcmb_GA = ft_appendtimelock(cfg, betas_plcmb{allsubs});
betas_bl_plcmb_GA = ft_appendtimelock(cfg, betas_bl_plcmb{allsubs});
clear betas_plcmb betas_bl_plcmb


%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
% neighbours           = ft_prepare_neighbours(cfg, tstat1_plCmb_GA); % define neighbouring channels
neighbours           = ft_prepare_neighbours(cfg, betas_plcmb_GA); % define neighbouring channels

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

stat = ft_timelockstatistics(cfg, betas_plcmb_GA, betas_bl_plcmb_GA);


% save
filename = sprintf('/project/3011085.02/results/stat_glm_gamma_time_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat', 'betas_plcmb_GA', 'betas_bl_plcmb_GA');

ft_diary('off')



