function erf_osc_analysis_erf(subj, isPilot)
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
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata_no_ica.mat', subj), 'dataClean');
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

% baseline correct with last 100ms from baseline window
cfg=[];
cfg.baseline = [-0.1 0];
dataOnsetLBl = ft_timelockbaseline(cfg, dataL);
dataOnsetMBl = ft_timelockbaseline(cfg, dataM);
dataOnsetRBl = ft_timelockbaseline(cfg, dataR);

cfg=[];
cfg.offset = -(dataL.trialinfo(:,5)-dataL.trialinfo(:,4)); % shift data with #sampleShift - #sampleOnset
dataShiftL = ft_redefinetrial(cfg, dataL); % shift uncorrected data to shift onset
dataShiftLBl = ft_redefinetrial(cfg, dataOnsetLBl); % shift corrected data to shift onset
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
dataShiftM = ft_redefinetrial(cfg, dataM);
dataShiftMBl = ft_redefinetrial(cfg, dataOnsetMBl);
cfg.offset = -(dataR.trialinfo(:,5)-dataR.trialinfo(:,4));
dataShiftR = ft_redefinetrial(cfg, dataR);
dataShiftRBl = ft_redefinetrial(cfg, dataOnsetRBl);

% baseline correct shifted data with 100ms preShift window
cfg=[];
cfg.baseline = [-0.1 0];
dataShiftLBlPre = ft_timelockbaseline(cfg, dataShiftL);
dataShiftMBlPre = ft_timelockbaseline(cfg, dataShiftM);
dataShiftRBlPre = ft_timelockbaseline(cfg, dataShiftR);



%% Time-lock analysis
% based on grating onset (baseline window corrected)
cfg=[];
cfg.vartrllength = 2;
cfg.channel = 'MEG';
% cfg.keeptrials='yes';
tlOnsetL = ft_timelockanalysis(cfg, dataOnsetLBl);
tlOnsetM = ft_timelockanalysis(cfg, dataOnsetMBl);
tlOnsetR = ft_timelockanalysis(cfg, dataOnsetRBl);
% based on grating shift (baseline window corrected)
tlShiftL = ft_timelockanalysis(cfg, dataShiftLBl);
tlShiftM = ft_timelockanalysis(cfg, dataShiftMBl);
tlShiftR = ft_timelockanalysis(cfg, dataShiftRBl);
% based on grating shift (100ms preShift corrected)
tlShiftPreL = ft_timelockanalysis(cfg, dataShiftLBlPre);
tlShiftPreM = ft_timelockanalysis(cfg, dataShiftMBlPre);
tlShiftPreR = ft_timelockanalysis(cfg, dataShiftRBlPre);

% downsample for clarity in plots
cfg=[];
cfg.resamplefs = 200;
tlOnsetL_rs = ft_resampledata(cfg, tlOnsetL);
tlOnsetM_rs = ft_resampledata(cfg, tlOnsetM);
tlOnsetR_rs = ft_resampledata(cfg, tlOnsetR);
tlShiftL_rs = ft_resampledata(cfg, tlShiftL);
tlShiftM_rs = ft_resampledata(cfg, tlShiftM);
tlShiftR_rs = ft_resampledata(cfg, tlShiftR);
tlShiftPreL_rs = ft_resampledata(cfg, tlShiftPreL);
tlShiftPreM_rs = ft_resampledata(cfg, tlShiftPreM);
tlShiftPreR_rs = ft_resampledata(cfg, tlShiftPreR);


%% save

save(sprintf('timelock_subj%d_no_ica', subj), 'tlOnsetL_rs', 'tlOnsetM_rs', 'tlOnsetR_rs',...
    'tlShiftL_rs', 'tlShiftM_rs', 'tlShiftR_rs', 'tlShiftPreL_rs', 'tlShiftPreM_rs', 'tlShiftPreR_rs')
%% Plot topography
%{
cfg=[];
cfg.layout = 'CTF275_helmet.mat';
cfg.latency = [0 0.2];
cfg.xlim=[0 0.2];
cfg.zlim=[-1e-13 1e-13];
cfg.colormap = 'jet';

% based on grating onset (baseline window corrected)
figure;
subplot(1,3,1);ft_topoplotER(cfg, tlOnsetL_rs);title('Left, onset (baseline window corrected)');
subplot(1,3,2);ft_topoplotER(cfg, tlOnsetM_rs);title('Middle, onset (baseline window corrected)');
subplot(1,3,3);ft_topoplotER(cfg, tlOnsetR_rs);title('Right, onset (baseline window corrected)');
% based on grating shift (baseline window corrected)
figure;
subplot(1,3,1);ft_topoplotER(cfg, tlShiftL_rs);title('Left, shift (baseline window corrected)');
subplot(1,3,2);ft_topoplotER(cfg, tlShiftM_rs);title('Middle, shift (baseline window corrected)');
subplot(1,3,3);ft_topoplotER(cfg, tlShiftR_rs);title('Right, shift (baseline window corrected)');
% based on grating shift (100ms preShift corrected)
figure;
subplot(1,3,1);ft_topoplotER(cfg, tlShiftPreL_rs);title('Left, shift (baseline 100ms pre)');
subplot(1,3,2);ft_topoplotER(cfg, tlShiftPreM_rs);title('Middle, shift (baseline 100ms pre)');
subplot(1,3,3);ft_topoplotER(cfg, tlShiftPreR_rs);title('Right, shift (baseline 100ms pre)');

%% plot timecourses

cfg=[];
cfg.xlim = [-0.5 1];
cfg.ylim = [-5e-13 3e-13];
figure;
% based on grating onset (baseline window corrected)
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlOnsetL_rs);title('Left, onset (baseline window corrected)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlOnsetM_rs);title('Middle, onset (baseline window corrected)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlOnsetR_rs);title('Right, onset (baseline window corrected)');
% based on grating shift (baseline window corrected)
figure;
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlShiftL_rs);title('Left, shift (baseline window corrected)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlShiftM_rs);title('Middle, shift (baseline window corrected)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlShiftR_rs);title('Right, shift (baseline window corrected)');
% based on grating shift (200ms preShift corrected)
figure;
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlShiftPreL_rs);title('Left, shift (baseline 100ms pre)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlShiftPreM_rs);title('Middle, shift (baseline 100ms pre)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlShiftPreR_rs);title('Right, shift (baseline 100ms pre)');
%}

end
