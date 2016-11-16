function erf_osc_freqstat

workSpace = whos;
diary('tmpDiary') % save command window output
for i = 1:numel(workSpace) % list all workspace variables
    eval(workSpace(i).name)
end

%% frequency statistics

cfg=[];
cfg.frequency = 'all';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistisc = 'maxsum';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;
design = [ones(1,nTrials), 2*ones(1,nTrials); 1:nTrials, 1:nTrials];
cfg.design = design;
cfg.ivar = 1;
cfg.uvar = 2;
cfg_neighb.method = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg_neighb, tfaHannLeft_onset);
[stat_on_bl] = ft_freqstatistics(cfg, tfaHannLeft_onset, tfaHannLeft_bl);
[stat_shift_bl] = ft_freqstatistics(cfg, tfaHannLeft_shift, tfaHannLeft_bl);

%% save
filename = sprintf('/project/3011085.02/Results/ERF_oscillation/freq/freqstat_subj%d', subj);
save(fullfile([filename '.mat']), 'stat_on_bl', 'stat_shift_bl'); 
diary off
movefile('tmpDiary', fullfile([filename, '.txt']))