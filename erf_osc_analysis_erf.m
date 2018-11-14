function erf_osc_analysis_erf(subj, erfoi)
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
if nargin<2 || isempty(erfoi)
    erfoi = 'reversal';
end


% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;
data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');

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

% realign to grating change or behavioral response
if strcmp(erfoi, 'reversal')
    cfg=[];
    cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
    data=ft_redefinetrial(cfg, data);
elseif strcmp(erfoi, 'motor')
    cfg=[];
    cfg.offset = -(data.trialinfo(:,6)-data.trialinfo(:,4));
    data=ft_redefinetrial(cfg, data);
end

%% Time-lock analysis
%select channel
cfg              = [];
cfg.vartrllength = 2;
cfg.channel      = {'MEG'};
cfg.keeptrials   = 'yes';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfilttype = 'firws';
cfg.preproc.lpfreq = 30;
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
tlck = ft_timelockanalysis(cfg, data);

cfg=[];
cfg.latency = [-0.1 0.5];
tlck = ft_selectdata(cfg, tlck);

%{
%planar gradient transformation
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tlck);

cfg.planarmethod    = 'sincos';
tl_planar      = ft_megplanar(cfg, tlck);

tl_plcmb       = ft_combineplanar([], tl_planar);

load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaPow');
[val, idx] = sort(gammaPow, 'descend');
qSize = round(length(gammaPow)/4);
cfg=[];
cfg.trials = idx(1:qSize);
cfg.avgoverrpt = 'yes';
if strcmp(erfoi, 'reversal')
    cfg.latency = [0 0.5];
elseif strcmp(erfoi, 'motor')
    cfg.latency = [-0.5 0];
end
q1 = ft_selectdata(cfg, tl_plcmb);
cfg.trials = idx(end-qSize+1:end);
q4 = ft_selectdata(cfg, tl_plcmb);
%}

%% save

    filename = sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_timelock_%s', subj,subj, erfoi);

save(fullfile([filename '.mat']), 'tlck', '-v7.3')
ft_diary('off')

end
