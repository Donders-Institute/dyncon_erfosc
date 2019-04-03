function erfosc_analysis_gamma_pow(subj, isPilot)
% This function makes an estimate of gamma power increase active-baseline
% in order to exclude participants that do not show a clear gamma peak

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
erfosc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
end
data=data.dataClean;
fs = data.fsample;

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

for iTrl=1:size(data.trial,2)
    blonset(iTrl,1) = data.time{iTrl}(1);
end
cfg=[];
% baseline window
% cfg.trl = [data.trialinfo(:,3), data.trialinfo(:,4) data.trialinfo(:,3)-data.trialinfo(:,4)];
% the first column represents the start of the baseline period. Its sample
% number is inaccurate w.r.t. the time axis in the data (possibly because
% of previous use of ft_selectdata?). This results in NaNs in the data. 
% Thus, don't use the samplenumber provided by trialinfo, but calculate on 
% the spot based on time axis.
cfg.trl = [data.trialinfo(:,4)+blonset*fs, data.trialinfo(:,4)];
cfg.trl = [cfg.trl, cfg.trl(:,1)-cfg.trl(:,2)];
dataBl = ft_redefinetrial(cfg, data);

cfg=[];
% post erf period till stimulus reversal
cfg.trl = [data.trialinfo(:,4)+0.4*fs data.trialinfo(:,5) 0.4*fs*ones(nTrials,1)];
dataAct = ft_redefinetrial(cfg, data);

cfg         = [];
cfg.length  = 0.5; % 0.5 second windows
cfg.overlap = 0.5; % half overlap
bl = ft_redefinetrial(cfg, dataBl); % PROBLEM!!
% bl.trial contains NaNs. in ft_freqanalysis, this results in the whole
% array becoming NaN. specifically in ft_preproc_polyremoval of
% ft_specest_mtmfft
act = ft_redefinetrial(cfg, dataAct);


%% gamma power
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.taper       = 'hanning';
cfg.foilim      = [2 100];
cfg.keeptrials  = 'no'; % average baseline over trials
cfg.channel     = {'MLO', 'MZO', 'MRO'};
powBl           = ft_freqanalysis(cfg, bl);
powAct          = ft_freqanalysis(cfg, act);

cfg=[];
cfg.parameter = 'powspctrm';
cfg.operation = 'x1./x2';
powRatio      = ft_math(cfg, powAct, powBl);
% cfg.operation = 'subtract';
% powDiff       = ft_math(cfg, powAct, powBl);
gamRange = [find(powRatio.freq==30) find(powRatio.freq==90)];
gammaAvg       = mean(powRatio.powspctrm,1);
[maxP maxIdx]  = max(gammaAvg(gamRange(1):gamRange(2)));
freqs = powRatio.freq(gamRange(1):gamRange(2));
peakFreq_gamma = freqs(maxIdx);

cfg=[];
cfg.frequency = [30 90];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';
gamRatio        = ft_selectdata(cfg, powRatio);
gamRatio        = gamRatio.powspctrm;

%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/sub-%03d_pow', subj, subj);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/sub-%03d_pow', subj, subj);
end
save(fullfile([filename '.mat']), 'powRatio', 'gamRatio', 'peakFreq_gamma');
ft_diary('off')

end

