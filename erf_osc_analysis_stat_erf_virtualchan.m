function erf_osc_analysis_stat_erf_virtualchan(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

cfg=[];
cfg.keeptrials = 'yes';
for subj=allsubs
    tmp1{subj} = load(sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_erf_virtualchan_%s.mat', subj, subj, 'reversal'),'tlck');
    tlck{subj} = ft_timelockanalysis(cfg, tmp1{subj}.tlck);
	tmp2{subj} = load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel.mat', subj, subj),'gammaPow'); % load gamma power
    g{subj} = tmp2{subj}.gammaPow;
end
clear tmp1 tmp2

%% Bin based on gamma power
for subj=allsubs
[~, idx{subj}] = sort(g{subj}, 'descend');
end
cfg=[];
cfg.avgoverrpt = 'yes';
for subj=allsubs
cfg.trials = idx{subj}(1:length(g{subj})/4);
tl_q1{subj} = ft_selectdata(cfg, tlck{subj});
cfg.trials = idx{subj}(end-length(g{subj})/4+1:end);
tl_q4{subj} = ft_selectdata(cfg, tlck{subj});
end


%% flip ambiguous dipole and splitt up trials
cfg=[];
cfg.avgoverrpt = 'yes';
for subj=allsubs
    tlck{subj} = ft_selectdata(cfg, tlck{subj});
end


% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
tlck_GA = ft_appendtimelock(cfg, tlck{allsubs});
q1_GA = ft_appendtimelock(cfg, tl_q1{allsubs});
q4_GA = ft_appendtimelock(cfg, tl_q4{allsubs});
clear tlck tl_q1 tl_q4

% dipole flip
erf = squeeze(tlck_GA.trial);
erf = diag(1./std(erf,[],2))*erf; % normalize before SVD

[u,s,v]=svd(erf);

u2 = repmat(sign(u(:,1)), [1, 1, length(tlck_GA.time)]);
q1_GA.trial = u2.*q1_GA.trial;
q4_GA.trial = u2.*q4_GA.trial;

%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.channel          = 'gam_pow';
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
filename = sprintf('/project/3011085.02/analysis/stat_erf_virtualchan_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat');

ft_diary('off')



