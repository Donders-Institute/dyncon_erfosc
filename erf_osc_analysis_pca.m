% 
% trialinfo columns:
% 1: trialnumber
% 2: position (-1=left, 0=middle, 1=right)
% 3: sample of baseline onset
% 4: sample of grating onset
% 5: sample of grating shift (=0 if no shift)
% 6: sample of response (=0 if no response or if response too early)
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
% ft_timelockbaseline gives NaNs to the baselinewindow. These give errors
% in PCA.
cfg=[];
cfg.baseline = [-0.5 0];
cfg.parameter = 'trial';
cfg.keeptrials = 'yes';
dataOnsetLBl = ft_timelockbaseline(cfg, dataL);
dataOnsetLBl.trial(isnan(dataOnsetLBl.trial))=0;
dataOnsetMBl = ft_timelockbaseline(cfg, dataM);
dataOnsetMBl.trial(isnan(dataOnsetMBl.trial))=0;
dataOnsetRBl = ft_timelockbaseline(cfg, dataR);
dataOnsetRBl.trial(isnan(dataOnsetRBl.trial))=0;



cfg=[];
cfg.offset = -(dataL.trialinfo(:,5)-dataL.trialinfo(:,4)); % shift data with #sampleShift - #sampleOnset
dataShiftLBl = ft_redefinetrial(cfg, dataOnsetLBl); % shift corrected data to shift onset
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
dataShiftMBl = ft_redefinetrial(cfg, dataOnsetMBl);
cfg.offset = -(dataR.trialinfo(:,5)-dataR.trialinfo(:,4));
dataShiftRBl = ft_redefinetrial(cfg, dataOnsetRBl);


cfg=[];
cfg.resamplefs = 200;
cfg.detrend = 'no';
cfg.demean = 'yes';
dataL_shiftRs = ft_resampledata(cfg, dataShiftLBl);
dataM_shiftRs = ft_resampledata(cfg, dataShiftMBl);
dataR_shiftRs = ft_resampledata(cfg, dataShiftRBl);

cfg=[];
cfg.channel = {'MLO', 'MZO', 'MRO', 'MLP', 'MZP', 'MRP'};
cfg.latency = [0 0.3];
dataL_shiftRs = ft_selectdata(cfg, dataL_shiftRs);
dataM_shiftRs = ft_selectdata(cfg, dataM_shiftRs);
dataR_shiftRs = ft_selectdata(cfg, dataR_shiftRs);

%% perform the independent component analysis (i.e., decompose the data)
cfg=[];
cfg.method          = 'pca';
% cfg.channel         = {'MLO', 'MZO', 'MRO', 'MLP', 'MZP', 'MRP'};
cfg.numcomponent    = 5;
[compM] = ft_componentanalysis(cfg, dataM_shiftRs);
% compM.unmixing contains the eigenvectors in the rows for the n largest
% eigenvalues.

%% Compute components scores
% how much of the variation in the data is explained by the principal
% components? If there are a few components that explain a lot of the
% variability, take these components en compute the inner product of the
% component with the single-trial waveform. This gives components scores,
% which can be used to divide the data into groups.
for iTrl = 1:nTrials
    
    for iComp=1:5
for iSample=1:size(dataM_shiftRs.trial,3)
    tmp(iSample,iComp) = dot(compM.unmixing(iComp,:), dataM_shiftRs.trial(iTrl,:,iSample)');
end
    end
scoreUnnormalized = sum(tmp,1);

S(iTrl,:) = 1e12*(1/size(dataM_shiftRs.trial,3))* scoreUnnormalized;
end
