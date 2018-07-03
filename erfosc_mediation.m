
erf_osc_datainfo;
cd /project/3011085.02/scripts/erfosc/analysis_JM_data
d = dir('*corrpowlcmv_mvepeaks3.mat');

k=1;
for subj = allsubs
    tmp = load(d(k).name,'pow', 'X');% load gamma power and ERF amplitude for single trials
    X{k} = tmp.X;
    pow{k} = tmp.pow;
    tmp = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj)); % load reaction times
    rt{k} = tmp.rt;
    k=k+1;
end

% Add mediation toolbox
addpath(genpath('/project/3011085.02/scripts/MediationToolbox'));
addpath(genpath('/project/3011085.02/scripts/CanlabCore'));
%% single parcel
% take, per subject, the ERF amplitude in the parcel with the highest 
% effect size (gamma-erf), of only those parcels contributing to this
% effect
%{
load('/project/3011085.02/results/stat_peakpicking3.mat')
y= S.rho(:,stat.mask);
[~,maxeffect_idx] = max(y');

chans = stat.label(stat.mask);
chanidx = match_str(stat.label, chans);

for k=1:32
    [result,stat1,stat2] = mediation(pow{k}, rt{k}, X{k}(chanidx(maxeffect_idx(k)),:)');
    ab(1,k) = stat1.t(5);
    cp(1,k) = stat1.t(3);
end

% how to test this?
% ttest(ab)?
%}

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
        [result, stat1,stat2] = mediation(pow{m}, rt{m}, X{m}(k,:)');
        ab(m,k) = stat1.t(5); % can be 'stat1.mean(5)' for beta weight instead of tval.
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

% QUESTION:
% I don't know if we can test this similarly as the correlation between
% gamma power and ERF amplitude?
cfgs = [];
cfgs.method='montecarlo';
cfgs.design=[ones(1,n) ones(1,n)*2;1:n 1:n];
cfgs.statistic='ft_statfun_wilcoxon';
cfgs.numrandomization=10000;
cfgs.ivar=1;
cfgs.uvar=2;
cfgs.parameter='ab';
cfgs.correctm='cluster';
cfgs.clusterthreshold='nonparametric_individual';
cfgs.connectivity = parcellation2connectivity(atlas);
cfgs.neighbours = cfgs.connectivity;
cfgs.correcttail = 'prob';
cd /project/3011085.02/lastround
stat=ft_freqstatistics(cfgs, source, ref);


%% plot
stat.brainordinate = atlas;
load cortex_inflated_shifted; atlas.pos=ctx.pos;
cfgx = [];
cfgx.method='surface';
cfgx.funparameter='stat';
cfgx.funcolormap = flipud(brewermap(64, 'RdBu'));
cfgx.funcolorlim = 'maxabs';
cfgx.maskstyle='colormix';
cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
stat.brainordinate=atlas;
ft_sourceplot(cfgx,stat)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([91 18])

