function erf_osc_analysis_glm_gamma_time(subj, isPilot)
% do a linear regression of pre-change gamma power over time.
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
%     gammaPow2(i) = (gammaChan.trial(i).pow);
end
gammaPow = gammaPow-mean(gammaPow);
nTrials = length(data_dss.trial);

[~, idxMax] = sort(gammaPow, 2, 'descend');

%% GLM on all trials
design = [ones(size(gammaPow)); gammaPow];
data=data_dss;
cfg=[];
cfg.lpifilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 40;
data = ft_preprocessing(cfg, data);
data.trial = cat(3,data.trial{:});
data.trial = permute(data.trial, [3,1,2]);
data.time = data.time{1};

if zscore
    t1 = nearest(data.time, (-0.5+1/fs));
    t2 = nearest(data.time, 0);
    S = std(data.trial(:,:,t1:t2), [], 3);
    S = repmat(S,[1,1,length(data.time)]);
    M = mean(data.trial(:,:,t1:t2),3);
    M = repmat(M,[1,1,length(data.time)]);
    data.trial=(data.trial-M)./S;
end

for k=1:length(data.label)
    Y = squeeze(data.trial(:,k,:));
    betas_trials(:,:,k) = design'\Y;
end


%% planar gradiant transformation of beta weights
% put beta weights in timelock structure
tl=[];
tl.avg    = squeeze(betas_trials(2,:,:))';
tl.time   = data_dss.time{1};
tl.dimord = 'chan_time';
tl.label  = data_dss.label;
tl.grad   = data_dss.grad;

% planar combination
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
cfg.planarmethod    = 'sincos';
tlPlanar            = ft_megplanar(cfg, tl);
cfg                 = [];
tlPlanarCmb         = ft_combineplanar(cfg,tlPlanar);


%% Save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_gamma_time', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_gamma_time', subj);
end
save(fullfile([filename '.mat']), 'betas_trials', 'tlPlanarCmb', '-v7.3');
diary off
movefile(diaryname, fullfile([filename '.txt']));

