function [data, data_shift] = erf_osc_preprocessing_eyedata(subj)
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
screenRight=1919;
screenLeft=0;
screenTop=0;
screenBottom=1079;

% range of voltages and recording range as defined in FINAL.INI
minVoltage=-5;
maxVoltage=5;
minRange=0;
maxRange=1;

totDist=780; % distance eye-screen in mm

% if you want to calculate the visual degree per pixel, physical dimensions
% of the pixels have to be known (so projected screen size).
screenWidth_y=270;
screenWidth_x=480; % has to be 'projected width' in mm (not actual width)

% convert ADC-voltref back to GAZE position
xGaze=[];
yGaze=[];
for iTrial=1:size(datatmp.trial,2);% these are trials with right-left voluntary saccades
    
    
    idxChanHor=find(strcmp(datatmp.label,'UADC005'));
    idxChanVer=find(strcmp(datatmp.label,'UADC006'));
    
    voltageHor=datatmp.trial{iTrial}(idxChanHor,:);
    voltageVer=datatmp.trial{iTrial}(idxChanVer,:);
    
    % see eyelink1000_UserManual for formulas:
    RangeHor = (voltageHor-minVoltage)./(maxVoltage-minVoltage); % voltage range proportion
    S_h = RangeHor.*(maxRange-minRange)+minRange; % proportion of screen width or height
    
    R_v = (voltageVer-minVoltage)./(maxVoltage-minVoltage);
    S_v = R_v.*(maxRange-minRange)+minRange;
    
    xGaze{iTrial} = S_h.*(screenRight-screenLeft+1)+screenLeft;
    yGaze{iTrial} = S_v.*(screenBottom-screenTop+1)+screenTop;
    
    
    % gaze position on center of screen (where fixation dot is)
    Xgaze_origin=xGaze;
    Ygaze_origin=yGaze;
    
    xGaze{iTrial} = xGaze{iTrial} - round(screenRight/2);%%presentation/PTB assumes 1 more pixel than there actually is.
    %presenation x-position = screen position + 1 pixel to the right:
    %presentation x-pos(0) = 1919/2
    
    yGaze{iTrial} = yGaze{iTrial} - round(screenBottom/2);
end
% Xgaze and Ygaze now are the eye positions on the screen (relative to
% fixation point). We can now convert them to visual angles (relative to
% fixation point).

% if you want to know the visual angle per pixel, calculate this
%
% X
screenResX=screenRight+1;
visAngX_rad = 2 * atan(screenWidth_x/2/totDist);
visAngX_deg = visAngX_rad * (180/pi);
visAngPerPixX = visAngX_deg / screenResX;

dotSize = 20; % fixation dot is 20 pixels wide (see concentric_grating_experiment)
visDegFixDot = visAngPerPixX * dotSize;

% Y
screenResY=screenBottom+1;
visAngY_rad = 2 * atan(screenWidth_y/2/totDist);
visAngY_deg = visAngY_rad * (180/pi);
visAngPerPixY = visAngY_deg / screenResY;
%}

% visual angle in reference to middle of the screen
for iTrial = 1:size(datatmp.trial,2)
    tmp1.trial{iTrial} = rad2deg(atan((xGaze{iTrial}/2)/(totDist*screenResX/screenWidth_x)))*2;
    tmp2.trial{iTrial} = rad2deg(atan((yGaze{iTrial}/2)/(totDist*screenResY/screenWidth_y)))*2;
end
data = ft_appenddata([], data, tmp1, tmp2);


%% load clean data
load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subjk, subj));
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

data = rmfield(data,'cfg');
data_shift = rmfield(data_shift, 'cfg');


save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_eyedata.mat', subj, subj),'data', 'data_shift');