function erf_osc_analysis_gamma_pow(subj, isPilot)
% 
% trialinfo columns:
% 1: trialnumber
% 2: position (-1=left, 0=middle, 1=right)
% 3: sample of baseline onset
% 4: sample of grating onset
% 5: sample of grating shift (=0 if no shift)
% 6: sample of response (=0 if no response or if response too early)
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


%% FFT, powerspectrum
cfg=[];
cfg.channel = {'MZO', 'MZP', 'MLO', 'MLP', 'MRO', 'MRP'};
cfg.latency = [-1+1/fs 0];
dataBaseline = ft_selectdata(cfg, data);
% take second preceding shift (NOTE: is it confounding if this includes
% grating onset, which has higher gamma peak freq?)
dataActive  = ft_selectdata(cfg, dataShift);

cfg=[];
cfg.foi = peakFreq;
cfg.method = 'mtmfft';
cfg.output = 'pow';
% cfg.tapsmofrq = 1;
cfg.taper = 'hanning';
cfg.keeptrials = 'yes';
cfg.pad = 1;
powActive = ft_freqanalysis(cfg, dataActive);
powBaseline = ft_freqanalysis(cfg, dataBaseline);

cfg=[];
cfg.operation = '(x1-x2)./x2';
cfg.parameter = 'powspctrm';
gammaPow = ft_math(cfg, powActive, powBaseline);

%% save
filename = sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_pow_%d', subj, subj);
save(fullfile([filename '.mat']), 'gammaPow');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));


end

