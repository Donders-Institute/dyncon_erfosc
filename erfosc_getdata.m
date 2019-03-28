function [data_onset, data_shift, data_response] = erfosc_getdata(data, comp, channel)
% Manipulates processed MEG data, including trial selection, time locking
% to different zero point and resampling.
%
% INPUT
%   data (struct): FieldTrip data structure
%   comp (struct): output of componentanalysis
%   channel (string): which channels to include in the selection
%       (default='MEG'). 
%
%   OUTPUT
%   data_onset: data, time locked to stimulus onset
%   data_shift: data, time locked to stimulus change
%   data_response: data, time locked to behavioral response

if nargin<3
  channel = 'MEG';
end

if nargin<2
  comp = [];
end

if ~isempty(comp)
  resamplefs = 150;
else
  resamplefs = 600;
end

data.cfg = []; % this clears a very large cfg, which can be recovered from file
fsorig   = data.fsample;

% select only shift trials, with valid response
idxM    = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,6)>data.trialinfo(:,5));
% idxM    = find(data.trialinfo(:,5)>0); % include trials independent of response, for neuroimage reviewer 2.
nTrials = length(idxM);

cfg         = [];
cfg.trials  = idxM;
cfg.channel = channel;
data        = ft_selectdata(cfg, data);
if ~isempty(comp)
  cfg  = rmfield(cfg, 'channel');
  comp = ft_selectdata(cfg, comp);
end

% find out which trials have response after end of trial, so you can
% exclude them
cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
data_reversal_tmp = ft_redefinetrial(cfg, data);
if ~isempty(comp)
  cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4)).*150./1200;
  comp = ft_redefinetrial(cfg, comp);
end

trlLatency = zeros(nTrials,1);
for iTrial=1:nTrials
    trlLatency(iTrial) = data_reversal_tmp.time{iTrial}(end);
end
idx_trials = find(trlLatency(:)>(data.trialinfo(:,6)-data.trialinfo(:,5))/1200);

cfg         = [];
cfg.trials  = idx_trials;
cfg.channel = channel;
data        = ft_selectdata(cfg, data);
if ~isempty(comp)
  cfg  = rmfield(cfg, 'channel');
  comp = ft_selectdata(cfg, comp);
end
clear data_reversal_tmp

% get the stim-onset aligned data
cfg         = [];
cfg.latency = [-0.75 0.75-1/fsorig];
data_onset  = ft_selectdata(cfg, data);

% get the response locked data
cfg           = [];
cfg.offset    = -(data.trialinfo(:,6)-data.trialinfo(:,4));
data_response = ft_redefinetrial(cfg, data);

% get the shift aligned data
cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
data       = ft_redefinetrial(cfg, data);

cfg         = [];
cfg.latency = [-0.75 0.75];
data_shift  = ft_selectdata(cfg, data);

if ~isempty(comp)
  ix = zeros(numel(comp.trial),2);
  for m = 1:numel(comp.trial)
    ix(m,1) = nearest(comp.time{m},-0.75);
    ix(m,2) = nearest(comp.time{m},0.5-1./150);
  end
  for m = 1:numel(comp.trial)
    comp.trial{m} = comp.trial{m}(:,ix(m,1):ix(m,2));
    comp.time{m}  = comp.time{m}(ix(m,1):ix(m,2));
  end
  cfg = [];
  cfg.time = comp.time;
  data_shift = ft_resampledata(cfg, data_shift);
  data_shift.trial = data_shift.trial + (comp.topo*comp.unmixing)*data_shift.trial;
  
  nsmp = cellfun('size',data_shift.trial, 2);
  ok   = nsmp==min(nsmp);
  sel  = find(~ok);
  for m = sel(:)'
    data_shift.trial{m} = data_shift.trial{m}(:,1:min(nsmp));
    data_shift.time{m}  = data_shift.time{m}(1:min(nsmp));
  end
  data_shift.time(1:end) = data_shift.time(1);
  
  data_onset = [];
  
else
  cfg         = [];
  cfg.resamplefs = resamplefs;
  data_onset  = ft_resampledata(cfg, data_onset);

  cfg         = [];
  cfg.resamplefs = resamplefs;
  data_shift  = ft_resampledata(cfg, data_shift);
end
