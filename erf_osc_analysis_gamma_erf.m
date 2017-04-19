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


% close all

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
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/dss.mat', subj), 'data_dss');
%     load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));    
end
fs=data_dss.fsample;
for i=1:length(gammaChan.trial)
    gammaPow(i) = log(gammaChan.trial(i).pow);
end
gammaPow = gammaPow-mean(gammaPow);
 nTrials = length(data_dss.trial);

 
%% GLM on binned trials
groupSize = round(nTrials/5);
[~, idxMax] = sort(gammaPow, 2, 'descend');

erfdata = data_dss;
% cfg=[];
% cfg.lpifilter = 'yes';
% cfg.lpfilterord = 2;
% cfg.lpfilttype = 'firws';
% erfdata = ft_preprocessing(cfg, erfdata);
erfdata.time=erfdata.time{1};
erfdata.trial=erfdata.trial(idxMax);
erfdata.trialinfo=erfdata.trialinfo(idxMax,:);
erfdata.trial = cat(3, erfdata.trial{:});
erfdata.trial = permute(erfdata.trial, [3,1,2]);


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
   cfg.trials = init:z;
   tmp = ft_selectdata(cfg, erfdata);
   erf{i} = ft_timelockanalysis(cfgb, tmp);
   var = erf{i}.var;
   cfg=[];
   cfg.baseline = [-0.05+1/fs 0];
   erf{i} = ft_timelockbaseline(cfg, erf{i});
   gamma(i) = mean(gammaPow(idxMax(init:z)));
   erf{i}.var = var;
   init=init+round(groupSize/2);
end

erfdata.trial = [];
for i = 1:round((nTrials/groupSize)*2)-1
erfdata.trial(i,:,:) = erf{i}.avg;
erfdata.var(i,:,:) = erf{i}.var;
end
% 
% t1 = nearest(erfdata.time, (-0.5+1/fs));
% t2 = nearest(erfdata.time, 0);
% S = sqrt(erfdata.var);
% M = mean(erfdata.trial(:,:,t1:t2),3);
% M = repmat(M,[1,1,length(erfdata.time)]);
% erfdata.trial=(erfdata.trial-M)./S;

design = [ones(size(gamma)); gamma];
for k=1:length(erfdata.label)
Y = squeeze(erfdata.trial(:,k,:));
betas_bin(:,:,k) = design'\Y;
end

%% GLM on all trials
design = [ones(size(gammaPow)); gammaPow];
data=data_dss;
data.trial = cat(3,data.trial{:});
data.trial = permute(data.trial, [3,1,2]);
data.time = data.time{1};

t1 = nearest(data.time, (-0.5+1/fs));
t2 = nearest(data.time, 0);
S = std(data.trial(:,:,t1:t2), [], 3);
S = repmat(S,[1,1,length(data.time)]);
M = mean(data.trial(:,:,t1:t2),3);
M = repmat(M,[1,1,length(data.time)]);
data.trial=(data.trial-M)./S;

for k=1:length(data.label)
Y = squeeze(data.trial(:,k,:));
betas_trials(:,:,k) = design'\Y;
end

%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/erf_gamma', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/erf_gamma', subj);
end
save(fullfile([filename '.mat']), 'betas_bin','betas_trials', 'erfdata', 'gamma');
diary off
movefile(diaryname, fullfile([filename '.txt']));

