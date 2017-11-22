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
    betas{subj} = rmfield(tmp{subj}.betas, 'elec');
    if strcmp(erfoi, 'motor')
        betas_bl{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time_%s.mat', subj, 'reversal'), 'betas_bl');
        betas_bl{subj} = rmfield(betas_bl{subj}.betas_bl, 'elec');
    else
        betas_bl{subj} = rmfield(tmp{subj}.betas_bl, 'elec');
    end
end
clear tmp

cfg=[];
cfg.method='template';
for subj=allsubs
    cfg.neighbours      = ft_prepare_neighbours(cfg, betas{subj});
    betas{subj} = ft_megplanar(cfg, betas{subj});
    betas_bl{subj} = ft_megplanar(cfg, betas_bl{subj});
end

for subj=allsubs
betas{subj} = ft_combineplanar([], betas{subj});
betas_bl{subj} = ft_combineplanar([], betas_bl{subj});
end
%% prepare for statistics

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
betas_GA = ft_appendtimelock(cfg, betas{allsubs});
betas_bl_GA = ft_appendtimelock(cfg, betas_bl{allsubs});
clear betas_plcmb betas_bl_plcmb

if strcmp(erfoi, 'motor')
    betas_bl_GA.time = betas_GA.time;
end

% take average baseline
betas_bl_GA.trial = repmat(mean(betas_bl_GA.trial,3), [1,1, length(betas_GA.time)]);
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

stat = ft_timelockstatistics(cfg, betas_GA, betas_bl_GA);
stat.cfg = rmfield(stat.cfg, 'previous');


% save
filename = sprintf('/project/3011085.02/results/stat_glm_gamma_time_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat');

ft_diary('off')



