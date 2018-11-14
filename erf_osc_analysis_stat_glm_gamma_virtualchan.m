function erf_osc_analysis_stat_glm_gamma_virtualchan(erfoi)

if nargin<1 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_glm_gamma_virtualchan_%s.mat', subj, subj, erfoi));
    betas{subj} = tmp{subj}.betas;
end
clear tmp

for subj=allsubs
    tmp{subj} = load(sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_erf_virtualchan_%s.mat', subj, subj, 'reversal'),'tlck');
    tlck{subj} = ft_timelockanalysis([], tmp{subj}.tlck);
end
clear tmp


%% prepare for statistics

% get grand average
cfg               = [];
cfg.appenddim     = 'rpt';
tlck_GA = ft_appendtimelock(cfg, tlck{allsubs});
betas_GA = ft_appendtimelock(cfg, betas{allsubs});
clear betas tlck

% dipole flip
erf = squeeze(tlck_GA.trial);
erf = diag(1./std(erf,[],2))*erf; % normalize before SVD

[u,s,v]=svd(erf);

u2 = repmat(sign(u(:,1)), [1, 1, length(tlck_GA.time)]);
betas_GA.trial = u2.*betas_GA.trial;
tlck_GA.trial = u2.*tlck_GA.trial;

% take zeros as reference; the beta weights are not biased.
ref = betas_GA;
ref.trial = ref.trial*0;

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

stat = ft_timelockstatistics(cfg, betas_GA, ref);
stat.cfg = rmfield(stat.cfg, 'previous');

% save
filename = sprintf('/project/3011085.02/analysis/stat_glm_gamma_virtualchan_%s', erfoi);
save(fullfile([filename, '.mat']), 'stat');

ft_diary('off')



