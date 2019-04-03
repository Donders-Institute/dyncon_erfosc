function erfosc_headlocation(subj)

if nargin<1
    subj = 4;
end
if isempty(subj)
    subj = 4;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
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

%% load data, define trials
erfosc_datainfo;
isPilot=true;

cfg=[];
if isPilot
    cfg.dataset = pilotsubjects(subj).dataset;
    cfg.logfile = load(pilotsubjects(subj).logfile);% load log file
else
    cfg.dataset = subjects(subj).dataset;
    cfg.logfile = load(subjects(subj).logfile);
end
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erfosc_trialfun';
cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);


%% calculate headposition
% preprocess the headposition data
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
    'HLC0021','HLC0022','HLC0023', ...
    'HLC0031','HLC0032','HLC0033'};
cfg.continuous = 'yes';

headpos = ft_preprocessing(cfg);

% calculate the mean coil position per trial
nTrials = length(headpos.sampleinfo);
for iTrl = 1:nTrials
    coil1(:,iTrl) = [mean(headpos.trial{1,iTrl}(1,:)); mean(headpos.trial{1,iTrl}(2,:)); mean(headpos.trial{1,iTrl}(3,:))];
    coil2(:,iTrl) = [mean(headpos.trial{1,iTrl}(4,:)); mean(headpos.trial{1,iTrl}(5,:)); mean(headpos.trial{1,iTrl}(6,:))];
    coil3(:,iTrl) = [mean(headpos.trial{1,iTrl}(7,:)); mean(headpos.trial{1,iTrl}(8,:)); mean(headpos.trial{1,iTrl}(9,:))];
end

cc = circumcenter(coil1, coil2, coil3);
cc_mean = mean(cc,2); % mean over trials
cc_dem = [cc - repmat(cc_mean,1,size(cc,2))]';% demean to obtain
% translations and rotations from the average position and orientation
cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';% get translations and
% rotations with reference to the first trial

meanposchange = mean(abs(cc_rel(:,1:3)*1000)) % from start-end trial

%% save
if isPilot
    filename = sprintf('/project/3011085.02/analysis/behavior/pilot-%03d/sub-%03d_pilot_headposition', subj, subj);
else
    filename = sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_headposition', subj, subj);
end
save(fullfile([filename '.mat']), 'cc', 'cc_dem', 'cc_rel');
diary off
movefile(diaryname, fullfile([filename, '.txt']));
end
