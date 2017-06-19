function erf_osc_analysis_stat_glm_gamma
% montecarlo based cluster statistics of the regression weights active vs
% baseline over all subjects.


% inititate diary
% ft_diary('on')

% load data
erf_osc_datainfo;

for subj=allsubs;
    w{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time.mat', subj),'tlPlanarCmb');
    w{subj} = w{subj}.tlPlanarCmb;
end

% select baseline
cfg1=[];
cfg1.avgovertime='no';
cfg1.latency = [-0.5 0];

% select active
cfg2=[];
cfg2.latency = [0 0.5];

cfg3=[];
cfg3.parameter = 'avg';
cfg3.operation = 'subtract';

for subj=allsubs
    bl{subj} = ft_selectdata(cfg1, w{subj});
    act{subj} = ft_selectdata(cfg2, w{subj});
    bl{subj}.time = act{subj}.time;
%     bl{subj}.avg = repmat(bl{subj}.avg, [1, length(bl{subj}.time)]);
    difference{subj} = ft_math(cfg3, act{subj}, bl{subj});
end
cfg=[];
cfg.keepindividual = 'yes';
actGA = ft_appendtimelock(cfg, act{allsubs});
blGA = ft_appendtimelock(cfg, bl{allsubs});

diffGA = ft_timelockgrandaverage([], difference{allsubs});


%% statistics
Nsub = length(allsubs);

cfg = [];
cfg.method      = 'template'; % try 'distance' as well
cfg.feedback    = 'no';
neighbours      = ft_prepare_neighbours(cfg, actGA); % define neighbouring channels

cfg = [];
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

% cfg.design
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, actGA, blGA);

% save
filename = '/project/3011085.02/results/stat_glm_gamma_time.mat';
save(filename, 'stat');

ft_diary('off')
