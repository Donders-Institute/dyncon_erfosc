function [data_onset, data_shift] = erfosc_preprocessing_eyedata(subj, megdata, dosave)
% Retrieve eyetracker data from the raw MEG data. From the eyedata those
% trials are selected that are present in the processed MEG data. X- and Y-
% gaze positions on the screen (in pixels) are transformed to X- and Y-
% positions in visual degrees, relative to the central fixation point.
% Eyedata is resampled at 600 Hz and time locked both to stimulus onset and
% stimulus change.
%
% INPUT
%   subj (int): subject ID, ranging from 1 to 33, excluding 10.
%
% OUTPUT
%   saves result on disk.
%   data (struct): eye-data containing 3 channels (x and y gaze, pupil
%       diameter (a.u.)). Time locked to stimulus onset.
%   data_shift (struct): eye-data containing 3 channels (x and y gaze, pupil
%       diameter (a.u.)). Time locked to stimulus change.

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(megdata)
  error('please give corresponding MEG data as input')
end
if nargin<3 || isempty(dosave)
    dosave = false;
end

%% Load data and define trials
erfosc_datainfo; % load subject specific info.

cfg=[];
cfg.dataset = subjects(subj).dataset;
load(subjects(subj).logfile);
cfg.logfile = subjects(subj).logfile;
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erfosc_trialfun';
cfg.trialdef.prestim = min(log.realBaselineDuration,log.setBaselineDuration);
cfg.trialdef.poststim = log.completeDurationGrating;
cfg.catchtrial = log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);

cfg.channel = {'UADC005', 'UADC006', 'UADC007'};
data_onset = ft_preprocessing(cfg);
datatmp = data_onset;


%% change X and Y positions into visual angle relative to fixation
tmp1  =rmfield(data_onset, {'trial', 'label'});
tmp1.label{1} = 'visAngleX';
tmp2=tmp1;
tmp2.label{1} = 'visAngleY';
[tmp1.trial, tmp2.trial] = transform_eyedata(datatmp);
cfg=[];
cfg.comment = 'use channels UADC005 and UADC006 to construct visAngleX: the distance from fixation in visual degrees.';
tmp1 = ft_annotate(cfg, tmp1);
cfg=[];
cfg.comment = 'use channels UADC005 and UADC006 to construct visAngleY: the distance from fixation in visual degrees.';
tmp2 = ft_annotate(cfg, tmp2);

data_onset = ft_appenddata([], data_onset, tmp1, tmp2);

%% load clean data
% get the correct trials
cfg=[];
cfg.trials = megdata.trialinfo(:,1);
data_onset = ft_selectdata(cfg, data_onset);

cfg=[];
cfg.resamplefs = 600;
data_onset = ft_resampledata(cfg, data_onset);

% get the shift aligned data
cfg        = [];
cfg.offset = -((data_onset.trialinfo(:,5)-data_onset.trialinfo(:,4)))/2;
data_shift = ft_redefinetrial(cfg, data_onset);

if dosave
    save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_eyedata.mat', subj, subj),'data', 'data_shift');
end