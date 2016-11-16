function erf_osc_analysis_gamma_freq(subj, isPilot)
% This function estimates the gamma peak frequency

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
    data = load(sprintf('/project/3011085.02/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;
fs   = data.fsample;

% select only shift trials, with valid response
idxM    = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM;
data       = ft_selectdata(cfg, data);


%% FFT, powerspectrum
cfg          = [];
cfg.channel  = {'MZO', 'MZP', 'MLO', 'MLP', 'MRO', 'MRP'};
cfg.latency  = [-1+1/fs 0];
dataBaseline = ft_selectdata(cfg, data);
cfg.latency  = [0.4+1/fs 1.75]; % take active time window after first erfs
dataActive   = ft_selectdata(cfg, data);

cfg            = [];
cfg.foi        = 30:1:100;
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.tapsmofrq  = 1;
cfg.taper      = 'dpss';
cfg.keeptrials = 'no';
cfg.pad        = 2;
powActive      = ft_freqanalysis(cfg, dataActive);
powBaseline    = ft_freqanalysis(cfg, dataBaseline);

cfg           = [];
cfg.operation = '(x1-x2)./x2';
cfg.parameter = 'powspctrm';
powDiff       = ft_math(cfg, powActive, powBaseline);

% average over channels, take the freq with max gamma pow diff
gammaAvg      = mean(powDiff.powspctrm,1);
[maxP maxIdx] = max(gammaAvg);
peakFreq      = powDiff.freq(maxIdx);


%% save
if ~exist(sprintf('/project/3011085.02/Results/ERF_oscillation/freq/%02d', subj), 'dir');
    mkdir(sprintf('/project/3011085.02/Results/ERF_oscillation/freq/%02d', subj));
end
filename = sprintf('/project/3011085.02/Results/ERF_oscillation/freq/%02d/gamma_peak_%d', subj, subj);
save(fullfile([filename '.mat']), 'peakFreq', 'gammaAvg');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));


end

