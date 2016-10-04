function erf_osc_analysis_gamma_phase(subj, isPilot)

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
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
load(sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_peak_%d', subj, subj), 'peakFreq');
data = data.dataClean;
fs = data.fsample;



% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
nTrials = length(idxM);

cfg=[];
cfg.trials = idxM(1:nTrials);
data = ft_selectdata(cfg, data);

cfg=[];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift = ft_redefinetrial(cfg, data);

cfg=[];
cfg.latency = [-0.5/peakFreq 0-1/fs];
dataPreShift = ft_selectdata(cfg, dataShift);

cfg=[];
cfg.resamplefs = 200;
data = ft_resampledata(cfg, data);

%% estimate gamma angle
cfg=[];
cfg.foi = peakFreq;
cfg.taper = 'hanning';
cfg.output = 'fourier';
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';
cfg.channel = {'MZO', 'MZP', 'MLO', 'MLP', 'MRO', 'MRP'};
fcomp = ft_freqanalysis(cfg, dataPreShift);
gammaAngle = angle(fcomp.fourierspctrm(:,1,:)); % in radians

%% save
filename = sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_angle_%d', subj, subj);
save(fullfile([filename '.mat']), 'gammaAngle');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));
end