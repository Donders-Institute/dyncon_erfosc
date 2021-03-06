function [tfa, baseline] = erf_osc_analysis_tfa(subj, isPilot, freqRange, zeropoint)
% Time frequency analysis.
%
% INPUT
%   subj (int): subject ID, ranging from 1 to 33, excluding 10 (default=1)
%   isPilot (logical): whether or not to apply on pilot data (default=0)
%   freRange (string): 'low' (default) or 'high', frequency range (below or 
%       above 30 Hz; affects configuration settings).
%   zeropoint (string): 'onset' (default) or 'reversal', what to time lock 
%       the data to. if 'onset' a seperate baseline estimate is saved.
%
% OUTPUT
%   saves result on disk.
%   tfa: time frequency estimate
%   baseline: baseline of time frequency estimate (only when
%       zeropoint='low').

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(freqRange)
    freqRange = 'low';
end
if nargin<4 || isempty(zeropoint)
    zeropoint = 'onset';% other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
end
data = data.dataClean;
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

%% Planar transformation
cfg                 = [];
cfg.method          = 'template';
cfg.template        = 'CTF275_neighb.mat';
cfg.neighbours      = ft_prepare_neighbours(cfg, data);
cfg.method          = 'sincos';
data_planar         = ft_megplanar(cfg, data);

%% TFA
cfg                  = [];
cfg.output           = 'pow';
cfg.channel          = 'MEG';
cfg.method           = 'mtmconvol';
cfg.toi              = -1.75:0.05:1.75;
cfg.keeptrials       = 'yes';

if strcmp(freqRange, 'low')
    cfg.taper        = 'hanning';
    cfg.pad          = 6;
    cfg.foi          = 2:2:30;% analysis 2 to 30 Hz in steps of 2 Hz
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    
elseif strcmp(freqRange, 'high')
    cfg.taper        = 'dpss';
    cfg.tapsmofrq    = 8; % 8 Hz freq smoothing on both sides
    cfg.foi          = 28:4:100;
    cfg.pad          = 6;
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*(1/4);
    
else
    error('freqRange must be *high* or *low*')
    keyboard
end

tfa_planar           = ft_freqanalysis(cfg,data_planar);

cfg                  = [];
cfg.method           = 'sum';
tfa                  = ft_combineplanar(cfg,tfa_planar);

%% get baseline
if strcmp(zeropoint, 'onset')
    cfg=[];
    cfg.latency = [-1 -0.25];
    baseline = ft_selectdata(cfg, tfa);
end

%% save

if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/sub-%03d_tfa_%s_%s', subj, subj, freqRange, zeropoint);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/sub-%03d_tfa_%s_%s', subj, subj, freqRange, zeropoint);
end
if strcmp(zeropoint, 'onset')
    save(fullfile([filename '.mat']), 'tfa', 'baseline');
elseif strcmp(zeropoint, 'reversal')
    save(fullfile([filename '.mat']), 'tfa');
end
ft_diary('off')

end

