function erf_osc_analysis_tfa(subj, isPilot, zeropoint)
%
% trialinfo columns:
% 1: trialnumber
% 2: position (-1=left, 0=middle, 1=right)
% 3: sample of baseline onset
% 4: sample of grating onset
% 5: sample of grating shift (=0 if no shift)
% 6: sample of response (=0 if no response or if response too early)
if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = false;
end
if isempty(isPilot);
    isPilot = false;
end
if nargin<3
    zeropoint = 'onset'; % other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if isempty(zeropoint);
    zeropoint = 'onset';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
end
data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
nTrials = length(idxM);

cfg=[];
cfg.trials = idxM;
cfg.channel = 'MEG';
data = ft_selectdata(cfg, data);

% data timelocked to grating shift
if strcmp(zeropoint, 'reversal')
    cfg=[];
    cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
    data = ft_redefinetrial(cfg, data);
end

% since trials are kept and are not the same length, memory wise it's more
% efficient to just take the length of the shortest trial.
cfg=[];
cfg.latency = [-1.75 1.75];
data = ft_selectdata(cfg, data);


%% TFA, low frequencies (2-30)

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 6;
cfg.foi          = 2:2:30;% analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -1.75:0.05:1.75;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials   = 'yes';
tfaLow = ft_freqanalysis(cfg, data);


%% TFA, high frequencies

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 8; % 8 Hz freq smoothing on both sides
cfg.foi          = 28:4:100;
cfg.pad          = 6;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*(1/4);
cfg.toi          = -1.75:0.05:1.75;
cfg.keeptrials   = 'yes';
tfaHigh = ft_freqanalysis(cfg,data);

%% get baseline
if strcmp(zeropoint, 'onset')
    cfg=[];
    cfg.latency = [-1 -0.25];
    cfg.avgoverrpt = 'yes';
    cfg.avgovertime = 'yes';
    baselineH = ft_selectdata(cfg, tfaHigh);
    
    cfg=[];
    cfg.latency = [-1 -0.25];
    cfg.avgoverrpt = 'yes';
    cfg.avgovertime = 'yes';
    baselineL = ft_selectdata(cfg, tfaLow);
end


%% save

if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/tfa_%s', subj, zeropoint);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s', subj, zeropoint);
end
if strcmp(zeropoint, 'onset')
    save(fullfile([filename '.mat']), 'tfaLow', 'tfaHigh', 'baselineL', 'baselineH');
elseif strcmp(zeropoint, 'reversal')
    save(fullfile([filename '.mat']), 'tfaLow', 'tfaHigh');
end
ft_diary('off')

end

