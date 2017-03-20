function erf_osc_analysis_corr(subj, isPilot)

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
% load gamma peak freq
if isPilot
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_peak.mat', subj), 'peakFreq');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_pow.mat', subj)); % load gamma power
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/aseo.mat', subj), 'amplitude', 'latency'); % load ERF
    % load(sprintf('/project/3011085.02/results/freq/pilot-0%d/gamma_angle_%d', subj, subj)); % load gamma phase
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_peak.mat', subj), 'peakFreq');
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_pow.mat', subj)); % load gamma power
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/aseo.mat', subj),  'amplitude', 'latency'); % load ERF
    % load(sprintf('/project/3011085.02/results/freq/subj-0%d/gamma_angle', subj)); % load gamma phase
end

peakFreq = 2*round(peakFreq/2);
gammaPow = gamPow.powspctrm;





%% gamma pow - ERF components
peakIdx = find(gamPow.freq==peakFreq);
[rAmp pAmp] = corr(gammaPow(:,peakIdx), amplitude, 'type', 'spearman')
[rLat pLat] = corr(gammaPow(:,peakIdx), latency, 'type', 'spearman')

% for iFreq = 1:length(gamPow.freq)
% [rA(iFreq,:), pA(iFreq,:)] = corr(gammaPow(:,iFreq), reconstructed.params.amplitude, 'type', 'spearman');
% [rL(iFreq,:), pL(iFreq,:)] = corr(gammaPow(:,iFreq), reconstructed.params.latency, 'type', 'spearman');
% end

%% gamma phase - ERF components
%{
gamAngle = rad2deg(gamAngle);
bins = 0:60:360;
angleBin = zeros(nTrials,1);
for iTrl = 1:nTrials
angleBin(iTrl) = bins(nearest(bins, gamAngle(iTrl)));
end
angleBin(angleBin==0) = 360;

idx60 = find(angleBin==60);
idx120 = find(angleBin==120);
idx180 = find(angleBin==180);
idx240 = find(angleBin==240);
idx300 = find(angleBin==300);
idx360 = find(angleBin==360);

gam60 = mean(gammaPow(idx60,1));
gam120 = mean(gammaPow(idx120,1));
gam180 = mean(gammaPow(idx180,1));
gam240 = mean(gammaPow(idx240,1));
gam300 = mean(gammaPow(idx300,1));
gam360 = mean(gammaPow(idx360,1));

gam = [gam60, gam120, gam180, gam240, gam300, gam360];

plot(60:60:360, gam, '.')
%}

%% save
% if isPilot
%     filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_pow', subj);
% else
%     filename = sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_pow', subj);
% end
% save(fullfile([filename '.mat']), 'gammaPow');
% diary off
% movefile(diaryname, fullfile([filename '.txt']));
end
