function erf_osc_analysis_tfa(subj, isPilot)
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
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response and the correct position
idxL = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==-1);
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
idxR = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==1);
nTrials = min([length(idxL) length(idxM) length(idxR)]);

cfg=[];
cfg.trials = idxL(1:nTrials); % equal number of trials for all conditions
dataL = ft_selectdata(cfg, data);
cfg.trials = idxM(1:nTrials);
dataM = ft_selectdata(cfg, data);
cfg.trials = idxR(1:nTrials);
dataR = ft_selectdata(cfg, data);

% data timelocked to grating shift
% cfg=[];
% cfg.offset = -(dataL.trialinfo(:,5)-dataL.trialinfo(:,4));
% dataShiftL = ft_redefinetrial(cfg, dataL);
% cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
% dataShiftM = ft_redefinetrial(cfg, dataM);
% cfg.offset = -(dataR.trialinfo(:,5)-dataR.trialinfo(:,4));
% dataShiftR = ft_redefinetrial(cfg, dataR);

% combine conditions
cfg=[];
cfg.trials = [idxL; idxM; idxR];
% cfg.latency = [-1 3.75]; % make all trials the same length
dataAll = ft_selectdata(cfg, data);

%% FFT, powerspectrum
cfg=[];
cfg.latency = [-0.5+1/fs 0];
dataBaseline = ft_selectdata(cfg, dataAll);
cfg.latency = [0.3+1/fs 0.8]; % take active time window after first erfs
dataActive  = ft_selectdata(cfg, dataAll);

cfg=[];
cfg.foilim = [2 100];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.tapsmofrq = 8;
cfg.taper = 'dpss';
cfg.channel = 'MEG';
cfg.keeptrials = 'yes';
powActive = ft_freqanalysis(cfg, dataActive);
powBaseline = ft_freqanalysis(cfg, dataBaseline);


%% TFA, low frequencies (2-30)

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 6;
cfg.foi          = 2:2:30;% analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -1:0.05:3.75;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.keeptrials   = 'yes';
tfaLow = ft_freqanalysis(cfg, dataAll);


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
cfg.toi          = -1.5:0.05:3.75;
cfg.keeptrials   = 'yes';
tfaHigh = ft_freqanalysis(cfg,dataAll);


save(sprintf('/home/electromag/matves/MATLAB/tfa_subj%d.mat', subj), 'powActive', 'powBaseline', 'tfaLow', 'tfaHigh', 'idxL', 'idxM', 'idxR');
%% frequency statistics
%{
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
[stat1] = ft_freqstatistics(cfg, tfaHannLeft_onset, tfaHannLeft_bl);
[stat2] = ft_freqstatistics(cfg, tfaHannLeft_shift, tfaHannLeft_bl);
%}

end

