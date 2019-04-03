function [data_onset, data_shift] = erfosc_preprocessing_eyedata(subj, dosave)
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
if nargin<2 || isempty(dosave)
    dosave = false;
end

%% Load data and define trials
erf_osc_datainfo; % load subject specific info.

cfg=[];
[~, name, ext] = fileparts(subjects(subj).dataset);
cfg.dataset = fullfile(['/project/3011085.02/erfosc/ses-meg01/', name, ext]);
cfg.logfile = load(sprintf('/project/3011085.02/erfosc/ses-beh01/sub%02dses01.mat', subj));
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erf_osc_mytrialfun';
cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);

cfg.channel = {'UADC005', 'UADC006', 'UADC007'};
data = ft_preprocessing(cfg);
datatmp = data;


%% change X and Y positions into visual angle relative to fixation
tmp1  =rmfield(data, {'trial', 'label'});
tmp1.label{1} = 'visAngleX';
tmp2=tmp1;
tmp2.label{1} = 'visAngleY';
[tmp1.trial, tmp2.trial] = transform_eyedata(datatmp);

data = ft_appenddata([], data, tmp1, tmp2);

%% load clean data
load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj));
[data_onset, ~, ~] = erfosc_getdata(dataClean, []);

% get the correct trials
cfg=[];
cfg.trials = data_onset.trialinfo(:,1);
data = ft_selectdata(cfg, data);

cfg=[];
cfg.resamplefs = 600;
data = ft_resampledata(cfg, data);

% get the shift aligned data
cfg        = [];
cfg.offset = -((data.trialinfo(:,5)-data.trialinfo(:,4)))/2;
data_shift = ft_redefinetrial(cfg, data);

data_onset = rmfield(data,'cfg');
data_shift = rmfield(data_shift, 'cfg');

if dosave
    save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_eyedata.mat', subj, subj),'data', 'data_shift');
end