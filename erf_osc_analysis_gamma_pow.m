function erf_osc_analysis_gamma_pow(subj, isPilot)
% This function estimates gamma power at the virtual channel (where gamma
% power was highest at gamma peak frequency)

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
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak.mat', subj), 'peakFreq');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gamPowData');
else
    load(sprintf('/project/3011085.02/results/freq/subj-%03d/gamma_peak.mat', subj), 'peakFreq');
    load(sprintf('/project/3011085.02/results/freq/subj-%03d/gamma_virtual_channel.mat', subj), 'gamPowData');
end
fs = gamPowData.fsample;

cfg                = [];
cfg.offset         = -(gamPowData.trialinfo(:,5)-gamPowData.trialinfo(:,4));
gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);

cfg          = [];
cfg.latency  = [-0.5+1/fs 0];
dataBl      = ft_selectdata(cfg, gamPowData);
dataPre     = ft_selectdata(cfg, gamPowDataShift);

peakFreq = 2*round(peakFreq/2);
smoothing = 6;
%% gamma power
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.tapsmofrq   = smoothing;
cfg.foilim      = [(peakFreq - 6*smoothing) (peakFreq + 6*smoothing)];
cfg.keeptrials  = 'no'; % average baseline over trials
gamPowPre       = ft_freqanalysis(cfg, dataBl);
cfg.keeptrials  = 'yes';
gamPowPost      = ft_freqanalysis(cfg, dataPre);

gamPow = gamPowPost;
gamPowPre.powspctrm = repmat(gamPowPre.powspctrm, [size(gamPow.powspctrm,1), 1]);

gamPow.powspctrm = (squeeze(gamPow.powspctrm) - gamPowPre.powspctrm)./gamPowPre.powspctrm;


%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_pow', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/subj-%03d/gamma_pow', subj);
end
save(fullfile([filename '.mat']), 'gamPow');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));


end

