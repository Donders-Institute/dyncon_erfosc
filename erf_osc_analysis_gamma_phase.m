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
diaryname = sprintf('tmpDiary_%s', datestr(now, 'dd.mm.yyyy_HH:MM:SS'));
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
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gamPowData');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak', subj), 'peakFreq');
else
    load(sprintf('/project/3011085.02/results/freq/subj-0%d/gamma_virtual_channel.mat', subj), 'gamPowData');
    load(sprintf('/project/3011085.02/results/freq/subj-0%d/gamma_peak', subj), 'peakFreq');
end
fs = gamPowData.fsample;
nTrials = length(gamPowData.trial);

cfg                = [];
cfg.offset         = -(gamPowData.trialinfo(:,5)-gamPowData.trialinfo(:,4));
gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);

cfg          = [];
cfg.latency  = [-5/peakFreq 0-1/fs]; % 5 cycles of peak freq
dataPreShift = ft_selectdata(cfg, gamPowDataShift);

%% estimate gamma angle
cfg            = [];
cfg.taper      = 'hanning';
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 5;
cfg.foilim     = [peakFreq peakFreq];
fcomp          = ft_freqanalysis(cfg, dataPreShift);
gamAngle       = angle(fcomp.fourierspctrm(:,1,:)); % in radians
% gamAngle       = radtodeg(gamAngle)+180; % shift to 0-360 degree

% bins = 0:60:360;
% angleBin = zeros(nTrials,1);
% for iTrl = 1:nTrials
% angleBin(iTrl) = bins(nearest(bins, gamAngle(iTrl)));
% end
% angleBin(angleBin==0) = 360;

%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_angle', subj);
else
    filename = sprintf('/project/3011085.02/results/freq/subj-%03d/gamma_angle', subj);
end
save(fullfile([filename '.mat']), 'gamAngle');
diary off
movefile(diaryname, fullfile([filename '.txt']));
end