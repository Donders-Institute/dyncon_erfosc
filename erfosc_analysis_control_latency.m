%% do everything independently per subject, do permutation statistics
clear
erf_osc_datainfo;

% dum=[];
for k=1:32
s = hanning(100)';
% load gammaPow and reaction times

    load(sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_rt.mat', allsubs(k), allsubs(k)), 'rt');
    load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/gamma_virtual_channel.mat', allsubs(k), allsubs(k)), 'gammaPow');

    g=gammaPow';

% model the signal
fs = 1200;
dt = 1/fs;

%determine the amount of delay in samples for each response.
delay = round(rt/dt);

maxlength = 1264; %max delay over all subjects
for n=1:numel(rt)
dsignal(n,:) = [zeros(1,delay(n)), s, zeros(1, maxlength-delay(n))];
end

design = [g, ones(size(g))];
x = dsignal-repmat(mean(dsignal,1), [numel(rt),1]);

cfg=[];
cfg.glm.statistic = 'beta';
cfg.glm.standardise = false;
tmp = statfun_glm(cfg, x', design');
dum{k} = dsignal;
% dum = [dum; dsignal];
beta{k} = tmp.stat(:,1);
clear dsignal
end

betas=[];
betas.trial(1,:,:) = cat(2, beta{:})';
betas.trial = permute(betas.trial, [2,1,3]);
betas.label{1} = 'chan1';
betas.time = dt:dt:size(betas.trial,3)*dt;
betas.dimord = 'rpt_chan_time';

ref = betas;
ref.trial = ref.trial*0;

%%
cfgs=[];
cfgs.statistic        = 'ft_statfun_depsamplesT';
cfgs.design(1,1:64) = [ones(1,32) 2*ones(1,32)];
cfgs.design(2,1:64) = [1:32 1:32];
cfgs.ivar=1;
cfgs.uvar=2;
cfgs.alpha = 0.05;
cfgs.method = 'montecarlo';
cfgs.clusterstatistic = 'maxsum';
cfgs.correctm = 'cluster';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';
cfgs.neighbours = [];
cfgs.parameter = 'trial';
cfgs.clusteralpha = 0.05;
stat=ft_timelockstatistics(cfgs, betas, ref);

save('/project/3011085.02/analysis/stat_control_simulation.mat', 'stat')