function erf_osc_analysis_gamma_phase(subj, isPilot)
% trialinfo columns:
% 1: trialnumber
% 2: position (-1=left, 0=middle, 1=right)
% 3: sample of baseline onset
% 4: sample of grating onset
% 5: sample of grating shift (=0 if no shift)
% 6: sample of response (=0 if no response or if response too early)

% select different conditions, take same number of trials, all trials
% contain shift and correct response. Baseline correct and shift time axis
% such that shift latency is zero point.
% * to be implemented *
% estimate gamma power in 1s preceding shift. select channel with highest
% gamma power in occipito-parietal sensors. In ERFs, estimate delay between
% shift presentation and onset ERF. Calculate gamma phase at that
% latency (which freq when it is broadband?). Look at distribution of gamma
% angle. See whether gamma angle predicts ERF peak latency/amplitude or
% reaction time.
if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
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
for i = 1:numel(workSpace) % list all workspace variables
    eval(workSpace(i).name)
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
data = data.dataClean;

% select only shift trials, with valid response and the correct position
idxL = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==-1);
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
idxR = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==1);
nTrials = min([length(idxL) length(idxM) length(idxR)]);

cfg=[];
cfg.trials = idxL(1:nTrials); % equal number of trials for all conditions
dataL = ft_selectdata(cfg, data);
cfg.trials = idxM(1:nTrials);
dataM = ft_selectdata(cfg, data);
cfg.trials = idxR(1:nTrials);
dataR = ft_selectdata(cfg, data);

% baseline correct with last 500ms from baseline window
cfg=[];
cfg.baseline = [-0.5 0];
dataOnsetLBl = ft_timelockbaseline(cfg, dataL);
dataOnsetMBl = ft_timelockbaseline(cfg, dataM);
dataOnsetRBl = ft_timelockbaseline(cfg, dataR);

cfg=[];
cfg.offset = -(dataL.trialinfo(:,5)-dataL.trialinfo(:,4)); % shift data with #sampleShift - #sampleOnset
dataShiftLBl = ft_redefinetrial(cfg, dataOnsetLBl); % shift corrected data to shift onset
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
dataShiftMBl = ft_redefinetrial(cfg, dataOnsetMBl);
cfg.offset = -(dataR.trialinfo(:,5)-dataR.trialinfo(:,4));
dataShiftRBl = ft_redefinetrial(cfg, dataOnsetRBl);



%%
cfg=[];
cfg.foi = 28:4:100;
cfg.taper = 'hanning';
cfg.output = 'fourier';
cfg.method = 'mtmfft';
cfg.keeptrials = 'yes';
cfg.channel = 'MEG'; %MZO01
fcomp = ft_freqanalysis(cfg,dataShiftLBl);
gammaAngle = angle(fcomp.fourierspctrm(:,1,:)); % in radians


diary off