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
close all

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
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
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

cfg=[];
cfg.latency = [-4.95 0.7]; % this will ensure the betas have the same size for all subjects
data = ft_selectdata(cfg, data);

for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
end


%% mutual information / GLM
% rt = (data.trialinfo(:,6)-data.trialinfo(:,5))'; 
design = [ones(size(gammaPow)); gammaPow-mean(gammaPow)];
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.latency = [0 0.5];
ft_topoplotER(cfgp, tlShift);

% channel = ;%input('what is the maximum channel?')
% load(sprintf('/project/3011085.02/results/erf/sub-%03d/erf_gamma.mat', subj), 'channel')

% ch_idx = find(strcmp(channel, data.label));
for k=1:length(data.label)
Y = squeeze(data.trial(:,k,:));
betas(:,:,k) = design'\Y;
end


%% Sort and bin trials on gamma power
groupSize = round(nTrials/5);
[~, idxMax] = sort(gammaPow, 2, 'descend');
% idxMax = idxMax(1:groupSize);
% [~, idxMin] = sort(gammaPow, 2, 'ascend');
% idxMin = idxMin(1:groupSize);
data.trial = data.trial(idxMax,:,:);
% data.dimord = 'rpt_chan_time';
meangamma = mean(gammaPow);
% cfg=[];
% cfg.trials = idxMax;
% cfg.avgoverrpt = 'yes';
% cfg.latency = [-0.1 inf];
% maxGam = ft_selectdata(cfg, data);
% cfg.trials = idxMin;
% minGam = ft_selectdata(cfg, data);
cfgb=[];
cfgb.vartrllength = 2;
init = 1;
for i = 1:round((nTrials/groupSize)*2)-1
    if init+groupSize-1>nTrials
        z = nTrials;
    else
        z = init+groupSize-1;
    end
   cfg=[];
   cfg.trials = idxMax(init:z);
   tmp = ft_selectdata(cfg, data);
   erf{i} = ft_timelockanalysis(cfgb, tmp);
   cfg=[];
   cfg.baseline = [-0.05+1/fs 0];
   erf{i} = ft_timelockbaseline(cfg, erf{i});
   erf{i}.gamma = mean(gammaPow(idxMax(init:z))-meangamma);
   init=init+round(groupSize/2);
end

erfdata=data;
erfdata.trial=[];
for i=1:length(erf)
    erfdata.trial{i} = erf{i}.avg;
    erfdata.var{i} = erf{i}.var;
    gamma(1,i) = erf{i}.gamma;
end


%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/erf_gamma', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/erf_gamma', subj);
end
save(fullfile([filename '.mat']), 'betas', 'erfdata', 'gamma');
diary off
movefile(diaryname, fullfile([filename '.txt']));

