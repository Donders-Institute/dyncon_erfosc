function erf_osc_analysis_stat_erf(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi='reversal';
end

ft_diary('on')

erf_osc_datainfo;

for subj=allsubs
    if strcmp(erfoi, 'reversal')
        tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock', subj), 'tlShift', 'tlShift_plcmb');
        tl{subj} = tmp{subj}.tlShift;
        tl_plcmb{subj} = tmp{subj}.tlShift_plcmb;
    else
        tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock', subj), 'tl', 'tl_plcmb');
        tl{subj} = tmp{subj}.tl;
        tl_plcmb{subj} = tmp{subj}.tl_plcmb;
    end
end
clear tmp


cfg=[];
cfg.appenddim='rpt';
tl_GA = ft_appendtimelock(cfg, tl{allsubs});
tl_plcmb_GA = ft_appendtimelock(cfg, tl_plcmb{allsubs});

% seperate active from baseline;
cfg=[];
cfg.latency = [-0.5 0];
tl_bl_GA = ft_selectdata(cfg, tl_GA);
tl_bl_plcmb_GA = ft_selectdata(cfg, tl_plcmb_GA);
cfg.latency = [0 0.5];
tl_act_GA = ft_selectdata(cfg, tl_GA);
tl_act_plcmb_GA = ft_selectdata(cfg, tl_plcmb_GA);

cfg=[];
cfg.parameter = 'trial';
cfg.operation = '0*x1';
zero_distribution = ft_math(cfg, tl_act_plcmb_GA);
%% statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template'; 
cfg.feedback         = 'no';
% neighbours           = ft_prepare_neighbours(cfg, tstat1_plCmb_GA); % define neighbouring channels
neighbours           = ft_prepare_neighbours(cfg, tl_act_plcmb_GA); % define neighbouring channels

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

stat = ft_timelockstatistics(cfg, tl_act_plcmb_GA, zero_distribution);
stat.cfg = rmfield(stat.cfg, 'previous');


% save
filename = sprintf('/project/3011085.02/results/stat_erf_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat', '-v7.3');

ft_diary('off')


