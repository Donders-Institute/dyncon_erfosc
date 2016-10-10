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
load(sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_virtual_channel_%d.mat', subj, subj), 'gamPowData');
load(sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_peak_%d', subj, subj), 'peakFreq');

cfg                = [];
cfg.offset         = -(gam_pow_data.trialinfo(:,5)-gam_pow_data.trialinfo(:,4));
gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);

cfg          = [];
cfg.latency  = [-0.5/peakFreq 0-1/fs];
dataPreShift = ft_selectdata(cfg, gamPowDataShift);

%% estimate gamma angle
cfg            = [];
cfg.foi        = peakFreq;
cfg.taper      = 'hanning';
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.keeptrials = 'yes';
fcomp          = ft_freqanalysis(cfg, dataPreShift);
gammaAngle     = angle(fcomp.fourierspctrm(:,1,:)); % in radians

%% save
filename = sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_angle_%d', subj, subj);
save(fullfile([filename '.mat']), 'gammaAngle');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));
end