
erf_osc_datainfo;
cd /project/3011085.02/analysis/

k=1;
for subj = allsubs
    filename = fullfile([datadir, 'corr/', sprintf('sub-%03d/sub-%03d_corrpowlcmv_gamma.mat', subj, subj)]);
    tmp = load(filename,'pow', 'X');% load gamma power and ERF amplitude for single trials
    X{k} = tmp.X;
    for m = 1:size(X{k},1)
        [~,~,rnkX{k}(m,:)] = unique(X{k}(m,:)); % rank transform ERF amplitude
    end
    pow{k} = tmp.pow;
    [~,~,rnkpow{k}] = unique(pow{k}); % rank transform gamma power
    tmp = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/sub-%03d_rt.mat', subj, subj)); % load reaction times
    rt{k} = tmp.rt;
    [~,~,rnkrt{k}] = unique(rt{k}); % rank transform RT
    k=k+1;
end

% Add mediation toolbox
addpath(genpath('/project/3011085.02/scripts/MediationToolbox'));
addpath(genpath('/project/3011085.02/scripts/CanlabCore'));

%% all parcels

% QUESTION:
% Wikipedia on Sobel test: when the mediator is included in a regression 
% analysis model with the independent variable, the effect of the 
% independent variable is reduced and the effect of the mediator remains 
% significant.

% in mediation.m from Canlab, the following are tested.
% columns of paths:
% -------------------------------------------------------------------------
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)

% Given the above, I am not sure which path we have to test. Is it just ab,
% which is the the mediation effect of M on X-->Y? Also, should we take the
% stat1.mean or stat1.t? As I understand, the values in 'mean' represent
% the beta coefficients. So do we continue with betas or with t-values?

% calculate mediation effect, seperately for every subject-parcel
% combination. Save the first level t-values.
for m=1:32
    for k=1:370
        [result, stat1,stat2] = mediation(rnkpow{m}, rnkrt{m}, rnkX{m}(k,:)');
        ab(m,k) = stat1.mean(5); % product of beta values a (X->M) and b (M->Y)
    end
end

load atlas_subparc374_8k.mat

% restructure into fieldtrip structure
source=[];
source.brainordinate = atlas;
source.label = atlas.parcellationlabel;
source.dimord = 'rpt_chan_freq';
source.freq = 0;
source.ab = zeros(32, 374);
exclude_label = match_str(source.brainordinate.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
idx = 1:374;
idx(exclude_label)=[];
source.ab(:,idx) = ab;
source.ab(:,exclude_label) = nan;
    
ref=source;
ref.ab(:)=0;
n=32;

gamma_erf_stat = load('/project/3011085.02/results/stat_peakpicking3.mat');

% not entirely sure

cfgs = [];
cfgs.method='montecarlo';
cfgs.design=[ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.statistic='ft_statfun_depsamplesT';
cfgs.numrandomization=10000;
cfgs.alpha = 0.05;
cfgs.clusteralpha = 0.05;
cfgs.ivar=1;
cfgs.uvar=2;
cfgs.parameter='ab';
cfgs.correctm='cluster';
cfgs.connectivity = parcellation2connectivity_midline(atlas);
cfgs.neighbours = cfgs.connectivity;
cfgs.correcttail = 'prob';
cfgs.clustertail = 1;
cfgs.tail = 1;
stat=ft_freqstatistics(cfgs, source, ref);

save('/project/3011085.02/results/stat_mediation_erf.mat', 'stat', 'source','rnkrt', 'rnkpow', 'rnkX')
%% plot
load cortex_inflated_shifted; atlas.pos=ctx.pos;
stat.brainordinate = atlas;

cfgx = [];
cfgx.method='surface';
cfgx.funparameter='stat';
cfgx.funcolormap = flipud(brewermap(64, 'RdBu'));
cfgx.funcolorlim = 'maxabs';
cfgx.maskstyle='colormix';
% cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
stat.brainordinate=atlas;
ft_sourceplot(cfgx,stat)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([91 18])

