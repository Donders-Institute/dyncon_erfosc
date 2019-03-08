function erf_osc_avgovertrl(subj, freqRange, erfoi)

load(sprintf('/project/3011085.02/results/freq/sub-%03d/sub-%03d_tfa_%s_%s.mat', subj, subj, freqRange, erfoi));

cfg=[];
cfg.avgoverrpt='yes';
tfa = ft_selectdata(cfg, tfa);

if strcmp(erfoi, 'onset')
    baseline = ft_selectdata(cfg, baseline);
    save(sprintf('/project/3011085.02/results/freq/sub-%03d/sub-%03d_tfa_%s_%s.mat', subj, subj, freqRange, erfoi), 'tfa', 'baseline', '-v7.3');
else
    save(sprintf('/project/3011085.02/results/freq/sub-%03d/sub-%03d_tfa_%s_%s.mat', subj, subj, freqRange, erfoi), 'tfa', '-v7.3');
end

end