function erf_osc_freqstat

workSpace = whos;
diaryname = tempname(fullfile([getenv('HOME'), '/tmp']));
diary(diaryname) % save command window output
fname = mfilename('fullpath')
datetime

fid = fopen(fullfile([fname '.m']));
tline = fgets(fid); % returns first line of fid
while ischar(tline) % at the end of the script tline=-1
    disp(tline) % display tline
    tline = fgets(fid); % returns the next line of fid
end
fclose(fid);

for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
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
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/freqstat', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/subj-%03d/freqstat', subj);
end
save(fullfile([filename '.mat']), 'stat_on_bl', 'stat_shift_bl'); 
diary off
movefile(diaryname, fullfile([filename, '.txt']))