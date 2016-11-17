function erf_osc_analysis_erf(subj, isPilot)
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
diary('tmpDiary') % save command window output
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
    data = load(sprintf('/project/3011085.02/Data/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/Data/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);

nTrials = length(idxM);

cfg=[];
cfg.trials = idxM;
cfg.channel = 'MEG';
dataM = ft_selectdata(cfg, data);

% baseline correct with last 100ms from baseline window
cfg=[];
cfg.baseline = [-0.1 0];
dataOnsetBl = ft_timelockbaseline(cfg, dataM);

cfg=[];
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
% dataShiftM = ft_redefinetrial(cfg, dataM);
dataShiftMBl = ft_redefinetrial(cfg, dataOnsetBl);

% baseline correct shifted data with 100ms preShift window
% cfg=[];
% cfg.baseline = [-0.1 0];
% dataShiftMBlPre = ft_timelockbaseline(cfg, dataShiftM);


%% Time-lock analysis
cfg=[];
cfg.vartrllength = 2;
cfg.channel = {'MLP', 'MRP'};
% cfg.keeptrials='yes';
% based on grating onset (baseline window corrected)
% tlOnset = ft_timelockanalysis(cfg, dataOnsetMBl);
% based on grating shift (baseline window corrected)
tlShift = ft_timelockanalysis(cfg, dataShiftMBl);
% based on grating shift (100ms preShift corrected)
% tlShiftPre = ft_timelockanalysis(cfg, dataShiftMBlPre);

% downsample for clarity in plots
cfg=[];
cfg.resamplefs = 200;
% tlOnsetM_rs = ft_resampledata(cfg, tlOnset);
tlShiftM_rs = ft_resampledata(cfg, tlShift);
% tlShiftPreM_rs = ft_resampledata(cfg, tlShiftPre);


%% save
if ~exist(sprintf('/project/3011085.02/Results/erf/%02d', subj), 'dir');
    mkdir(sprintf('/project/3011085.02/Results/erf/%02d', subj));
end
filename = sprintf('/project/3011085.02/Results/erf/%02d/timelock_%d', subj, subj);
save(fullfile([filename '.mat']), 'tlShiftM_rs')
diary off
movefile('tmpDiary', fullfile([filename, '.txt']));

end
