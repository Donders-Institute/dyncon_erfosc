function erf_osc_analysis_gamma_erf(subj, isPilot)

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

%% Initiate Diary
workSpace = whos;
diaryname = sprintf('/project/3011085.02/scripts/erfosc/tmpDiary_%s.txt', datestr(now, 'dd.mm.yyyy_HH:MM:SS.FFF'));
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
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
end

data = data.dataClean;
cfg         = [];
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

fs = data.fsample;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM;
data       = ft_selectdata(cfg, data);

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift  = ft_redefinetrial(cfg, data);

cfg=[];
cfg.baseline = [-0.050+1/fs 0];
data = ft_timelockbaseline(cfg, dataShift);


for i=1:length(gammaChan.trial)
    gammaPow(i) = gammaChan.trial(i).pow;
end

%% Sort and bin trials on gamma power
groupSize = round(nTrials/4);
[~, idxMax] = sort(gammaPow, 2, 'descend');
idxMax = idxMax(1:groupSize);
[~, idxMin] = sort(gammaPow, 2, 'ascend');
idxMin = idxMin(1:groupSize);

cfg=[];
cfg.trials = idxMax;
maxGam = ft_selectdata(cfg, data);
cfg.trials = idxMin;
minGam = ft_selectdata(cfg, data);


%% Timelock analysis
cfg=[];
cfg.vartrllength = 2;
tlMax = ft_timelockanalysis(cfg, maxGam);
tlMin = ft_timelockanalysis(cfg, minGam);




%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/erf_gamma', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/erf_gamma', subj);
end
save(fullfile([filename '.mat']), 'tlMax', 'tlMin');
diary off
movefile(diaryname, fullfile([filename '.txt']));

