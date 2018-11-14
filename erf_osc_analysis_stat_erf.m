function erf_osc_analysis_stat_erf(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi='reversal';
end

ft_diary('on')

erf_osc_datainfo;

for subj=allsubs
        tmp{subj} = load(sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_timelock_%s.mat', subj,subj, erfoi), 'q1', 'q4');
        q1{subj} = tmp{subj}.q1;
        q4{subj} = tmp{subj}.q4;
end
clear tmp

cfg=[];
cfg.appenddim='rpt';
q1_GA = ft_appendtimelock(cfg, q1{allsubs});
q4_GA = ft_appendtimelock(cfg, q4{allsubs});


%% statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
neighbours           = ft_prepare_neighbours(cfg, q1_GA); % define neighbouring channels

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

stat = ft_timelockstatistics(cfg, q1_GA, q4_GA);
stat.cfg = rmfield(stat.cfg, 'previous');


% save
filename = sprintf('/project/3011085.02/analysis/stat_erf_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat', '-v7.3');

ft_diary('off')


