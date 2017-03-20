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
    subj = 3;
end
if isempty(subj)
    subj = 3;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

%% initiate diary
workSpace = whos;
diaryname = sprintf('/project/3011085.02/scripts/erfosc/tmpDiary_%s.txt', datestr(now, 'dd.mm.yyyy_HH:MM:SS'));
diary(diaryname) % save command window output
fname = mfilename('fullpath')
datetime

fid = fopen(fullfile([fname '.m']));
tline = fgets(fid); % returns first line of fid
while ischar(tline) % at the end of the script tline=-1
    disp(tline) % display tline
    tline = fgets(fid); % returns the next line of fid
end
fclose(fid);

for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
end

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
fs = data.fsample;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
nTrials = length(idxM);

cfg=[];
cfg.trials = idxM(1:nTrials);
data = ft_selectdata(cfg, data);

% data timelocked to grating shift
% cfg=[];
% cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
% dataShiftM = ft_redefinetrial(cfg, dataM);

%% FFT, powerspectrum
cfg=[];
cfg.latency = [-0.5+1/fs 0];
dataBaseline = ft_selectdata(cfg, data);
cfg.latency = [0.3+1/fs 0.8]; % take active time window after first erfs
dataActive  = ft_selectdata(cfg, data);

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
tfaLow = ft_freqanalysis(cfg, data);


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
tfaHigh = ft_freqanalysis(cfg,data);


%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/tfa', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/tfa', subj);
end
save(fullfile([filename '.mat']), 'powActive', 'powBaseline', 'tfaLow', 'tfaHigh');
diary off
movefile(diaryname, fullfile([filename '.txt']));


end

