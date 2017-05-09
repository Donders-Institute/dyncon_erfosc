function erf_osc_analysis_dss(subj, isPilot)

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
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
end

data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM;
cfg.channel = 'MEG';
data       = ft_selectdata(cfg, data);

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift  = ft_redefinetrial(cfg, data);

cfg=[];
cfg.latency = [-0.05 0.5];
dataShift2 = ft_selectdata(cfg, dataShift);
cfg.latency = [-2.5 0.65]; % this will ensure the betas have the same size for all subjects
dataShift = ft_selectdata(cfg, dataShift);


cfg=[];
cfg.baselinewindow = [-0.050+1/fs 0];
cfg.demean = 'yes';
dataShift = ft_preprocessing(cfg, dataShift);
dataShift2 = ft_preprocessing(cfg, dataShift2);

%% DSS component analysis

windowGuess = [0 0.4];
state.X = 1;
nComp = 25;
% run a dss decomposition
params      = [];
% params.time = dataShift.time;
params.time = dataShift2.time;

params.demean = 'prezero';
params.pre = 0.05*fs-1;
params.pst = fs*windowGuess(1,2);
% [~,~,avgorig] = denoise_avg2(params,dataShift.trial,state);
[~,~,avgorig] = denoise_avg2(params,dataShift2.trial,state);

cfg          = [];
cfg.method   = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params = params;
cfg.dss.wdim = 100;
cfg.numcomponent = nComp;
cfg.cellmode = 'yes';
% comp = ft_componentanalysis(cfg, dataShift);

comp = ft_componentanalysis(cfg, dataShift2);

%% select components

tl=ft_timelockanalysis([],comp); % compute trial-average
[u,s,v]=svd(cov(tl.avg')); % singular value decomposition of dss components
% s contains the singular values; make it into a vector:
s=diag(s);
s=s./s(1); % normalize singular values by largest value
nComp_keep = length(find(s>=0.05)); % keep components that explain at least 5% of the variance

%% Backproject components to channel level data
cfg=[];
cfg.unmixing = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_orig = ft_componentanalysis(cfg, dataShift);

cfg=[];
cfg.components = nComp_keep+1:nComp;
data_dss = ft_rejectcomponent(cfg, comp_orig, dataShift);


%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/dss', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/dss2', subj);
end
save(fullfile([filename '.mat']), 'data_dss', 'nComp_keep', '-v7.3');
diary off
movefile(diaryname, fullfile([filename '.txt']));

