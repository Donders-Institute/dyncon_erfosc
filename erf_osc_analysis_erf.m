function erf_osc_analysis_erf(subj, isPilot)
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
    isPilot = false;
end
if isempty(isPilot);
    isPilot = false;
end

%% initiate diary
workSpace = whos;
diaryname = tempname(fullfile([getenv('HOME'), '/tmp']));
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

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
fs = data.fsample;
nTrials = length(idxM);

cfg=[];
cfg.trials = idxM;
cfg.channel = 'MEG';
dataM = ft_selectdata(cfg, data);

% realign to grating change
cfg=[];
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
dataShift = ft_redefinetrial(cfg, dataM);

% baseline correct with last 50 ms just before change
cfg=[];
cfg.baseline = [-0.050+1/fs 0];
data = ft_timelockbaseline(cfg, dataShift);

%% Time-lock analysis
cfg=[];
cfg.vartrllength = 2;
cfg.channel = {'MEG'};
tlShift = ft_timelockanalysis(cfg, data);


%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/timelock', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/timelock', subj);
end
save(fullfile([filename '.mat']), 'tlShift')
diary off
movefile(diaryname, fullfile([filename, '.txt']));

end
