function erf_osc_analysis_alpha_freq(subj, isPilot)
% This function estimates the gamma peak frequency

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
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
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
cfg.foi        = 7:1:14;
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
alphaAvg       = mean(powDiff.powspctrm,1);
[maxP maxIdx]  = max(alphaAvg);
peakFreq_alpha = powDiff.freq(maxIdx);


%% save
if isPilot
    filename = sprintf('/project/3011085.02/analysis/freq/pilot-%03d/sub-%03d_alpha_peak', subj, subj);
else
    filename = sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_alpha_peak', subj, subj);
end
save(fullfile([filename '.mat']), 'peakFreq_alpha', 'alphaAvg');
diary off
movefile(diaryname, fullfile([filename '.txt']));


end

