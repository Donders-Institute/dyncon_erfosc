function erf_osc_analysis_gamma_pow(subj, isPilot)
% This function estimates gamma power at the virtual channel (where gamma
% power was highest at gamma peak frequency)

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
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak.mat', subj), 'peakFreq_gamma');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gamPowData');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_peak.mat', subj), 'peakFreq_gamma');
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gamPowData');
end
fs = gamPowData.fsample;

cfg                = [];
cfg.offset         = -(gamPowData.trialinfo(:,5)-gamPowData.trialinfo(:,4));
gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);

cfg          = [];
cfg.latency  = [-0.5+1/fs 0];
dataBl      = ft_selectdata(cfg, gamPowData); % baseline
dataChange     = ft_selectdata(cfg, gamPowDataShift); % pre-change

peakFreq_gamma = 2*round(peakFreq_gamma/2);
smoothing = 6;
%% gamma power
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.tapsmofrq   = smoothing;
cfg.foilim      = [(peakFreq_gamma - 6*smoothing) (peakFreq_gamma + 6*smoothing)];
cfg.keeptrials  = 'no'; % average baseline over trials
gamPowBl       = ft_freqanalysis(cfg, dataBl);
cfg.keeptrials  = 'yes';
gamPowChange      = ft_freqanalysis(cfg, dataChange);

gamPow = gamPowChange;
gamPowBl.powspctrm = repmat(gamPowBl.powspctrm, [size(gamPow.powspctrm,1), 1]);

gamPow.powspctrm = (squeeze(gamPow.powspctrm) - gamPowBl.powspctrm)./gamPowBl.powspctrm;


%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_pow', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_pow', subj);
end
save(fullfile([filename '.mat']), 'gamPow');
diary off
movefile(diaryname, fullfile([filename '.txt']));


end

