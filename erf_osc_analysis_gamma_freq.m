function erf_osc_analysis_gamma_freq(subj, isPilot)
% This function estimates the gamma peak frequency

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

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;
fs   = data.fsample;

% select only shift trials, with valid response
idxM    = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM;
data       = ft_selectdata(cfg, data);


%% FFT, powerspectrum
cfg          = [];
cfg.channel  = {'MZO', 'MZP', 'MLO', 'MLP', 'MRO', 'MRP'};
cfg.latency  = [-1+1/fs 0];
dataBaseline = ft_selectdata(cfg, data);
cfg.latency  = [0.4+1/fs 1.75]; % take active time window after first erfs
dataActive   = ft_selectdata(cfg, data);

cfg            = [];
cfg.foi        = 30:1:90;
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.tapsmofrq  = 1;
cfg.taper      = 'dpss';
cfg.keeptrials = 'no';
cfg.pad        = 2;
powActive      = ft_freqanalysis(cfg, dataActive);
powBaseline    = ft_freqanalysis(cfg, dataBaseline);

cfg           = [];
cfg.operation = 'x1./x2';
cfg.parameter = 'powspctrm';
powRatio       = ft_math(cfg, powActive, powBaseline);

% average over channels, take the freq with max gamma pow diff
gammaAvg       = mean(powRatio.powspctrm,1);
[maxP maxIdx]  = max(gammaAvg);
peakFreq_gamma = powRatio.freq(maxIdx);


%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_peak', subj);
end
save(fullfile([filename '.mat']), 'peakFreq_gamma', 'gammaAvg');
ft_diary('off')


end

