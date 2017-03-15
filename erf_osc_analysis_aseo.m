function erf_osc_analysis_aseo(subj, isPilot)
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
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/timelock.mat', subj));
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
idxM        = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
cfg         = [];
cfg.trials  = idxM;
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

% redifine zero point to shift onset
cfg         = [];
cfg.offset  = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift   = ft_redefinetrial(cfg, data);

cfg         = [];
cfg.latency = [-0.05+1/fs 0.75];
dataShift   = ft_selectdata(cfg, dataShift);

cfg=[];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.05+1/fs 0];
dataShift = ft_preprocessing(cfg, dataShift);

%% Select time window
cfg=[];
cfg.baseline = [-0.05+1/fs 0];
cfg.baselinetype = 'abs';
tl = ft_timelockbaseline(cfg, tlShift);

t1 = nearest(tl.time, (-1+1/fs));
t2 = nearest(tl.time, 0);
s = std(tl.avg(:,t1:t2), [], 2);
m = mean(tl.avg(:, t1:t2), 2);
M = repmat(m, [1, size(tl.avg,2)]);

% two z scoring options: x/s or (x-m)/s?
y=tl;
y.avg=diag(1./s)*(tl.avg-M);
z = ft_globalmeanfield([], y);

y2=tl;
y2.avg = diag(1./s)*(tl.avg);
z2 = ft_globalmeanfield([], y2);

figure; subplot(2,1,1); plot(z.time, z.avg); title('(X-M)/s'); xlim([-0.05 0.5]);
subplot(2,1,2); cfg=[]; cfg.xlim = [-0.05 0.5]; ft_singleplotER(cfg, tl);
window = input('What is the latency window of every individual peak?');
% window = [0.056 0.094; 0.094 0.1517; 0.155 0.2117];
for i=1:size(window, 1)
    tmp = window(i,2)-window(i,1);
    tmp = tmp/4;
    windowvar(i,:) = [-tmp tmp];
end

%% Select channels
% repeat this for all peaks
for iPeak = 1:size(window, 1)
    
    % select critical window 0:300 ms
    cfg=[];
    cfg.latency = [window(iPeak,1) window(iPeak,2)];
    cfg.channel = {'MRP', 'MLP', 'MZP', 'MZO', 'MLO', 'MRO'};
    critWin = ft_selectdata(cfg, tlShift);
    
    % take mean of every channel
    avgChan = mean(critWin.avg, 2);
    nchan = length(critWin.label);
    
    % find 8 channels with max abs amplitude in critical window, then seperate
    % them into negative and positive (pole). Average them accordingly after
    % 'flipping' the negative pole
    selChan = 8;
    [~, maxAbsIdx] = sort(abs(avgChan), 1, 'descend');
    maxAbsIdx = maxAbsIdx(1:selChan);
    chanIdxPos = maxAbsIdx(find(avgChan(maxAbsIdx)>0));
    chanIdxNeg = maxAbsIdx(find(avgChan(maxAbsIdx)<0));
    
    cfg=[];
    cfg.channel = critWin.label(chanIdxPos);
    dataPos = ft_selectdata(cfg, dataShift);
    cfg.channel = critWin.label(chanIdxNeg);
    dataNeg = ft_selectdata(cfg, dataShift);
    %
    cfg = [];
    cfg.operation = '-x1';
    cfg.parameter = 'trial';
    tmp = ft_math(cfg, dataNeg);
    
    erfdata = ft_appenddata([], dataPos, tmp);
    cfg=[];
    cfg.avgoverchan = 'yes';
    erfdata = ft_selectdata(cfg, erfdata);
    
    cfg=[];
    cfg.xlim = [-0.05 0.5];
    figure; ft_singleplotER(cfg, erfdata);
    
    %% Analysis of Singletrial ERP and Ongoing acivity
    % Seperately model the evoked components and the ongoing activity from
    % a first estimate. Iterate this n times.
    
    cfg                      = [];
    cfg.method               = 'aseo';
    cfg.aseo.waveformInitSet = window
    cfg.aseo.numiteration    = 1;
    cfg.aseo.jitter          = windowvar
    [reconstructed{iPeak}, residual{iPeak}] = ft_singletrialanalysis(cfg, erfdata);
end
for iPeak = 1:size(window,1)
    latency(:, iPeak) = reconstructed{iPeak}.params.latency(:, iPeak);
    amplitude(:, iPeak) = reconstructed{iPeak}.params.amplitude(:, iPeak);
end

%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/aseo', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/aseo', subj);
end
save(fullfile([filename '.mat']), 'erfdata', 'reconstructed', 'residual', 'latency', 'amplitude')
diary off
movefile(diaryname, fullfile([filename, '.txt']));


end
