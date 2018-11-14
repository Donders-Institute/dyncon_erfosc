function erf_osc_resample(subj)

load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj));

cfg=[];
cfg.resamplefs = 600;
cfg.detrend = 'no';
dataClean = ft_resampledata(cfg, dataClean);

save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean', 'nTrialsPostClean','nTrialsPreClean', 'nTrialsValid', '-v7.3');
end
