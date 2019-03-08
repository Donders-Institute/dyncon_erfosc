function erf_osc_analysis_eyemovement(subj)

erf_osc_datainfo; % load subject specific info.
load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
cfg=[];
cfg.dataset = subjects(subj).dataset;%sprintf('/project/3011085.02/raw/sub-%03d/ses-meg01/sub11ses01_3011085.02_20170109_01.ds', subj);%subjects(subj).dataset;
cfg.logfile = load(sprintf('/project/3011085.02/raw/sub-%03d/ses-beh01/sub%02dses01.mat', subj, subj));
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erf_osc_mytrialfun';
cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
prestim = cfg.trialdef.prestim;
cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
poststim = cfg.trialdef.poststim;
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);
data = ft_preprocessing(cfg);

trlidx = find(ismember(data.sampleinfo(:,1), dataClean.sampleinfo(:,1)));
offset = prestim(trlidx);
stimulus = poststim(trlidx);

cfg=[];
cfg.channel = {'UADC005', 'UADC006'};
data_eye = ft_selectdata(cfg, data);
data_eye = rmfield(data_eye, 'cfg');
save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_eyedata.mat', subj, subj), 'data_eye', 'offset', 'stimulus');

%% select trials and redefine time axis

idxM = find(dataClean.trialinfo(:,5)>0 & dataClean.trialinfo(:,6)>0 & dataClean.trialinfo(:,6)>dataClean.trialinfo(:,5));
nTrials = length(idxM);

cfg=[];
cfg.trials = idxM;
data_eye = ft_selectdata(cfg, data_eye);
cfg.channel = 'MEG';
dataClean = ft_selectdata(cfg, dataClean);

% find out which trials have response after end of trial, so you can
% exclude them
cfg=[];
cfg.offset = -(dataClean.trialinfo(:,5)-dataClean.trialinfo(:,4));
data_reversal_tmp = ft_redefinetrial(cfg, dataClean);

for iTrial=1:nTrials
    trlLatency(iTrial) = data_reversal_tmp.time{iTrial}(end);
end
idx_trials = find(trlLatency'>((dataClean.trialinfo(:,6)-dataClean.trialinfo(:,5))/1200));
idx_trials_invalid = find(trlLatency'<((dataClean.trialinfo(:,6)-dataClean.trialinfo(:,5))/1200));

cfg=[];
cfg.trials = idx_trials;
data_eye = ft_selectdata(cfg, data_eye);
cfg.channel = 'MEG';
dataClean = ft_selectdata(cfg, dataClean);
clear data_reversal_tmp

cfg        = [];
cfg.offset = -(dataClean.trialinfo(:,5)-dataClean.trialinfo(:,4)); % trialinfo is specified in 1200 Hz. If data is resampled, it has to be taken care of for ft_redefinetrial.
dataShift  = ft_redefinetrial(cfg, dataClean);
dataShift_eye = ft_redefinetrial(cfg, data_eye);

cfg=[];
cfg.channel = 'UADC005';
X = ft_selectdata(cfg, dataShift_eye);
cfg.channel = 'UADC006';
Y = ft_selectdata(cfg, dataShift_eye);

save(sprintf('/project/3011085.02/analysis/eye/sub%03d.mat', subj), 'X', 'Y');





















