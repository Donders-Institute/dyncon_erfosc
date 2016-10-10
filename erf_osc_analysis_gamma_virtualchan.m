function erf_osc_analysis_gamma_virtualchan(subj, isPilot, isSegmentedMri)
% This function estimates gamma power at the estimated gamma peak frequency

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
if nargin<3
    isSegmentedMri = 0;
end
if isempty(isSegmentedMri)
    isSegmentedMri = 0;
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
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
load(sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_peak_%d', subj, subj), 'peakFreq');
data = data.dataClean;
cfg         = [];
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0 & data.trialinfo(:,2)==0);
nTrials = length(idxM);

cfg        = [];
cfg.trials = idxM(1:nTrials);
data       = ft_selectdata(cfg, data);

cfg        = [];
cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift  = ft_redefinetrial(cfg, data);

% cfg            = [];
% cfg.resamplefs = 200;
% data           = ft_resampledata(cfg, data);
% dataShift      = ft_resampledata(cfg, dataShift);
fs = data.fsample;


cfg         = [];
cfg.latency = [-1+1/fs 0];
dataPre     = ft_selectdata(cfg, data);
% take second preceding shift (NOTE: is it confounding if this includes
% grating onset, which has higher gamma peak freq?)
dataPost    = ft_selectdata(cfg, dataShift);
dataAll     = ft_appenddata([], dataPost, dataPre);

%% Frequency analysis
% calculate power in complete data (for calculation of common filter)
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 5;
cfg.foilim    = [peakFreq peakFreq];
freqAll       = ft_freqanalysis(cfg, dataAll);

% calculate power pre and post stimulus
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 5;
cfg.foilim    = [peakFreq peakFreq];
freqPre       = ft_freqanalysis(cfg, dataPre);
freqPost      = ft_freqanalysis(cfg, dataPost);

%% Source analysis
if isSegmentedMri
    load(pilotsubjects(subj).segmentedmri);
else
    mri = ft_read_mri(pilotsubjects(subj).mri);
    polhemus = ft_read_headshape(pilotsubjects(subj).headshape);
    
    % align mri with MEG data through polhemus headshape
    cfg                         = [];
    cfg.coordsys                = 'ctf';
    cfg.parameter               = 'anatomy';
    cfg.viewresult              = 'yes';
    cfg.method                  = 'headshape';
    cfg.headshape.headshape     = polhemus;
    cfg.headshape.interactive   = 'no';%'yes';
    cfg.headshape.icp           = 'yes';
    mri                         = ft_volumerealign(cfg, mri);
    
    % segment mri
    cfg             = [];
    cfg.write       = 'no';
    cfg.coordsys    = 'neuromag';
    [segmentedmri]  = ft_volumesegment(cfg, mri);
    save(pilotsubjects(subj).segmentedmri, 'segmentedmri');
end

% constructs a volume conduction model from the geometry of the head.
cfg         = [];
cfg.method  = 'singleshell';
headmodel   = ft_prepare_headmodel(cfg, segmentedmri);
headmodel   = ft_convert_units(headmodel, 'cm');

% constructs a source model, for example a 3-D grid or a cortical sheet.
cfg                 = [];
cfg.grid.resolution = 0.5;
cfg.grid.unit       = 'cm';
cfg.headmodel       = headmodel;
cfg.grad            = data.grad;
sourcemodel         = ft_prepare_sourcemodel(cfg);

% computes the forward model for many dipole locations on a regular 2D or
% 3D grid and stores it for efficient inverse modelling
cfg                 = [];
cfg.grad            = data.grad;
cfg.headmodel       = headmodel;
cfg.reducerank      = 2;
cfg.grid.resolution = 0.5;   % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
sourcemodel_lf      = ft_prepare_leadfield(cfg);

% calculate common filter
cfg                   = [];
cfg.grad              = freqAll.grad;
cfg.method            = 'dics';
cfg.frequency         = peakFreq;
cfg.grid              = sourcemodel_lf;
cfg.headmodel         = headmodel;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '5%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';
sourceAll             = ft_sourceanalysis(cfg, freqAll);

cfg.grid.filter       = sourceAll.avg.filter;
sourcePre             = ft_sourceanalysis(cfg, freqPre);
sourcePost            = ft_sourceanalysis(cfg, freqPost);
sourceDiff            = sourcePost;
sourceDiff.avg.pow    = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;

%% Create virtual channel of maximal gamma power
% the following position contains the max gamma power difference
[maxval, maxpowindx] = max(sourceDiff.avg.pow);
sourceDiff.pos(maxpowindx, :)
% we will create a virtual channel based on this location. In order to do
% this, we have to do timelockanalysis and use an LCMV beamformer, because
% this will pass the activity at the location of interest with unit gain,
% while optimally reducing activity from other locations. This filter then
% can be applied to the original MEG data

cfg                   = [];
cfg.covariance        = 'yes';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
tlock                 = ft_timelockanalysis(cfg, dataAll);

cfg                 = [];
cfg.method          = 'lcmv';
cfg.headmodel       = headmodel;
cfg.grid.pos        = sourcemodel.pos(maxpowindx, :);
cfg.grid.unit       = sourcemodel.unit;
cfg.lcmv.keepfilter = 'yes';
cfg.projectmom      = 'yes';
cfg.fixedori        = 'yes';
source_idx          = ft_sourceanalysis(cfg, tlock);

beamformerGamPow = source_idx.avg.filter;

gamPowData              = [];
gamPowData.label        = {'gam_pow'};
gamPowData.time         = data.time;
gamPowData.datainfo     = data.trialinfo;
for i=1:length(dataShift.trial)
    gamPowData.trial{i} = beamformerGamPow{1} * data.trial{i};
end

%% save
filename = sprintf('/home/electromag/matves/Results/ERF_oscillation/freq/%02d/gamma_virtual_channel_%d', subj, subj);
save(fullfile([filename '.mat']), 'gamPowData', 'beamformerGamPow');
diary off
movefile('tmpDiary', fullfile([filename '.txt']));


end

