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
load(sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_virtual_channel_%d.mat', subj, subj), 'gamPowData');

cfg                = [];
cfg.offset         = -(gam_pow_data.trialinfo(:,5)-gam_pow_data.trialinfo(:,4));
gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);

cfg          = [];
cfg.latency  = [-1+1/fs 0];
dataPre      = ft_selectdata(cfg, gamPowData);
dataPost     = ft_selectdata(cfg, gamPowDataShift);

%% gamma power
cfg             = [];
cfg.method      = 'mtmfft';
cfg.output      = 'pow';
cfg.tapsmofrq   = 5;
cfg.foilim      = [peakFreq peakFreq];
cfg.keeptrials  = 'yes';
gamPowPre       = ft_freqanalysis(cfg, dataPre);
gamPowPost      = ft_freqanalysis(cfg, dataPost);

gamPowDif           = gamPowPre;
gamPowDif.powspctrm = (gamPowPost.powspctrm-gamPowPre.powspctrm)./gamPowPre.powspctrm;


%% save
filename = sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_pow_%d', subj, subj);
save(fullfile([filename '.mat']), 'gamPowDif');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));


end

