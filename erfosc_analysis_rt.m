function [rt, idx_trials] = erfosc_analysis_rt(subj, dosave)
% calculate reaction times from trigger information and save on disk.
%
% INPUT
%   subj (int): subject ID, ranging from 1 to 33, excluding 10.
%
% OUTPUT
%   saves data on disk.
%   rt: reaction times (ms), after stimulus change
%   idx_trials: trial numbers of reaction time estimates

if nargin<1 || isempty(subj)
    subj=1;
end
if isempty(dosave)
    dosave=true;
end


%% load data, select trials, estimate RT
data=load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');

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

rt = (data.trialinfo(:,6)-data.trialinfo(:, 5))/data.fsample;
cfg=[];
cfg.comment = 'define rt by running "rt = (data.trialinfo(:,6)-data.trialinfo(:, 5))/data.fsample;".';
rt = ft_annotate(cfg, rt);

% save
filename = sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_rt', subj,subj);
if dosave
    save(fullfile([filename '.mat']), 'rt');
end


end