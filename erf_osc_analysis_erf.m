function erf_osc_analysis_erf(subj, isPilot)
% trialinfo columns:
% 1: trialnumber
% 2: position (-1=left, 0=middle, 1=right)
% 3: sample of baseline onset
% 4: sample of grating onset
% 5: sample of grating shift (=0 if no shift)
% 6: sample of response (=0 if no response or if response too early)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
end
data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
    idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,6)>data.trialinfo(:,5));
    nTrials = length(idxM);
    
    cfg=[];
    cfg.trials = idxM;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    
    % find out which trials have response after end of trial, so you can
    % exclude them
    cfg=[];
    cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
    data_reversal_tmp = ft_redefinetrial(cfg, data);
    
    for iTrial=1:nTrials
        trlLatency(iTrial) = data_reversal_tmp.time{iTrial}(end);
    end
    idx_trials = find(trlLatency'>((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
    idx_trials_invalid = find(trlLatency'<((data.trialinfo(:,6)-data.trialinfo(:,5))/1200));
    
    cfg=[];
    cfg.trials = idx_trials;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    clear data_reversal_tmp

% realign to grating change
cfg=[];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift = ft_redefinetrial(cfg, data);

% baseline correct with last 50 ms just before change
cfg=[];
cfg.baseline = [-0.050+1/fs 0];
data = ft_timelockbaseline(cfg, data);
dataShift = ft_timelockbaseline(cfg, dataShift);

%% Time-lock analysis
cfg              = [];
cfg.vartrllength = 2;
cfg.channel      = {'MEG'};
tl               = ft_timelockanalysis(cfg, data);
tlShift          = ft_timelockanalysis(cfg, dataShift);

% planar gradient transformation
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);

cfg.planarmethod    = 'sincos';
tl_planar           = ft_megplanar(cfg, tl);
tlShift_planar      = ft_megplanar(cfg, tlShift);

tl_plcmb            = ft_combineplanar([], tl_planar);
tlShift_plcmb       = ft_combineplanar([], tlShift_planar);



%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/timelock', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/timelock', subj);
end
save(fullfile([filename '.mat']), 'tl_plcmb', 'tlShift_plcmb', 'tl', 'tlShift')
ft_diary('off')

end
