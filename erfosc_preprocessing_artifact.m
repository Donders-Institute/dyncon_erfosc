function [dataClean, artfctdef, comp, headmotion, nTrials] = erfosc_preprocessing_artifact(subject, artfctdef, comp, visDegOffFixation)
% This function is interactive and can not be automated and/or
% paralellized.
% This function loads data from the ERF-oscillation experiment, segments
% it into trials and preprocesses (e.g. highpass filter, bandstop filter
% linenoise etc). Next, by visual inspection bad channels/trials are
% rejected and jump, muscle and EOG artifacts are selected. ICA components
% decomposition is aplied and components belonging to cardiac activity or
% eye artifacts are rejected from the data. The data is then saved on disk.
%
% use:
% subj, specify subject number (default = 1)
% isPilot, true if datasets belong to pilot subjects (default = false)
% existArtifact, true if artficats were already selected and saved (default
% = false)
% visDegOffFixation: threshold for deleting trials based on excessive eye
% movements (default: delete trials if deviation>5 visual degrees).
% global ft_default;
% ft_default = [];
% ft_default.reproducescript = datestr(now, 30);

if nargin<1 || isempty(subject)
  error('please specify subject as subject=subjects(x);')
end
if nargin<2 || isempty(artfctdef) || ~isstruct(artfctdef)
  existArtifact = false; else  existArtifact = true; end
if nargin<3 || isempty(comp) || ~isstruct(comp)
  doica = true; else doica=false; end
if nargin<4 || isempty(visDegOffFixation)
  visDegOffFixation = 5;
end


%% Load data and define trials
cfg=[];
cfg.dataset = subject.dataset;
cfg.logfile = subject.logfile;
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'erfosc_trialfun';
load(cfg.logfile);
cfg.trialdef.prestim = min(log.realBaselineDuration, log.setBaselineDuration);
cfg.trialdef.poststim = log.completeDurationGrating;
cfg.catchtrial = log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);
cfg2=cfg; % use this for finding EOG artifacts.
nTrials.preClean = size(cfg.trl,1);

% preprocess data
cfg.continuous = 'yes';
cfg.demean = 'yes';
cfg.padding = 8; % pad the data so the filter is the same for all trials (filter depends on trial length)
cfg.dftfilter = 'yes';
cfg.dftfreq   = [50+(-1:1)./cfg.padding 100+(-1:1)./cfg.padding 150+(-1:1)./cfg.padding 200+(-1:1)./cfg.padding]; % remove line noise
cfg.usefftfilt = 'yes'; % frequency multiplication of filter instead of convolution (slow)
cfg.hpfilter = 'yes';
cfg.hpfreq = 1; % remove slow drifts
cfg.hpfilttype = 'firws';
data = ft_preprocessing(cfg);



%%%%%%%%%%%%%%%%%%%%%%
% Artifact Selection %
%%%%%%%%%%%%%%%%%%%%%%
if ~existArtifact
  %% Reject bad channels/trials
  % open a GUI with summary info in order to reject bad channels and/or
  % trials with high variance.
  cfg          = [];
  cfg.method   = 'summary';
  cfg.channel = {'MEG'};
  cfg.layout = 'CTF275_helmet.mat'; % so timecourses and topographies of individual trials can be looked at
  dataRejVis        = ft_rejectvisual(cfg, data);
  % you can check the rejected trial number by typing
  artfctdef.badtrial=[];
  for i=1:size(dataRejVis.cfg.artfctdef.summary.artifact,1)
    artfctdef.badtrial(i)=find(data.sampleinfo(:,1)==dataRejVis.cfg.artfctdef.summary.artifact(i,1));
  end;
  artfctdef.badtrial = [artfctdef.badtrial, subject.badtrials];
  if artfctdef.badtrial
    disp(artfctdef.badtrial);
    data.badtrials = artfctdef.badtrial;
  else
    data.badtrials = [];
  end
  %% Find jump artifacts
  sprintf('find jump artifacts')
  cfg = [];
  % channel selection, cutoff and padding
  cfg.artfctdef.zvalue.channel    = 'MEG';
  cfg.artfctdef.zvalue.cutoff     = 20; % artifact threshold
  cfg.artfctdef.zvalue.trlpadding = 0;
  cfg.artfctdef.zvalue.artpadding = 0;
  cfg.artfctdef.zvalue.fltpadding = 0;
  % algorithmic parameters
  cfg.artfctdef.zvalue.cumulative    = 'yes';
  cfg.artfctdef.zvalue.medianfilter  = 'yes';
  cfg.artfctdef.zvalue.medianfiltord = 9;
  cfg.artfctdef.zvalue.absdiff       = 'yes';
  % make the process interactive
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  sprintf('find jump artifacts')
  [~, artfctdef.jump.artifact] = ft_artifact_zvalue(cfg, data);
  
  %% Find muscle artifacts
  cfg            = [];
  % channel selection, cutoff and padding
  cfg.artfctdef.zvalue.channel = 'MEG';
  cfg.artfctdef.zvalue.cutoff      = 4; % artifact threshold
  cfg.artfctdef.zvalue.trlpadding  = 0; % pad at both sides of trial till it is ... seconds.
  cfg.artfctdef.zvalue.fltpadding  = 0;
  cfg.artfctdef.zvalue.artpadding  = 0.2; % take 0.2s padding of an artifact (of the part that is above threshold)
  % algorithmic parameters
  cfg.artfctdef.zvalue.bpfilter    = 'yes';
  cfg.artfctdef.zvalue.bpfreq      = [110 140];
  cfg.artfctdef.zvalue.bpfiltord   = 9;
  cfg.artfctdef.zvalue.bpfilttype  = 'but';
  cfg.artfctdef.zvalue.hilbert     = 'yes';
  cfg.artfctdef.zvalue.boxcar      = 0.2;
  % make the process interactive
  cfg.artfctdef.zvalue.interactive = 'yes';
  
  sprintf('find muscle artifacts')
  [~, artfctdef.muscle.artifact] = ft_artifact_zvalue(cfg, data);
  
  %% find EOG blink artifacts
  % channel selection, cutoff and padding
  cfg2.artfctdef.zvalue.channel     = {'UADC006'}; % Eye tracker data
  cfg2.artfctdef.zvalue.cutoff      = 2;
  cfg2.artfctdef.zvalue.trlpadding  = 0;
  cfg2.artfctdef.zvalue.artpadding  = 0.01;
  cfg2.artfctdef.zvalue.fltpadding  = 0;
  % algorithmic parameters
  cfg2.artfctdef.zvalue.bpfilter   = 'yes';
  cfg2.artfctdef.zvalue.bpfilttype = 'but';
  cfg2.artfctdef.zvalue.bpfreq     = [1 30];
  cfg2.artfctdef.zvalue.bpfiltord  = 4;
  cfg2.artfctdef.zvalue.hilbert    = 'yes';
  % feedback
  cfg2.artfctdef.zvalue.interactive = 'yes';
  
  sprintf('find blink artifacts')
  [~, artfctdef.eyeblink.artifact] = ft_artifact_zvalue(cfg2);
  
  %% Find EOG saccade artifacts
  % this algorithm marks events as artifacts if eye movements are made
  % with more than 5 visual degree off the central fixation point.
  [artfctdef.eyepos.visAngleX, artfctdef.eyepos.visAngleY, visAngPerPixX, visAngPerPixY] = transform_eyedata(data);
  dotSize = 20; % fixation dot is 20 pixels wide (see concentric_grating_experiment)
  visDegFixDot = visAngPerPixX * dotSize;
  % NOTE: the rest of the eye movement artifact detection is done after
  % ICA based on visAngleX and Y.
end

%% Reject artifacts from complete data
cfg = [];
cfg.artfctdef = removefields(artfctdef, 'eyeblink2', 'eyesaccade'); % remove these after ICA
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-1 3.75];
dataNoArtfct = ft_rejectartifact(cfg, data);

% reject extra channels if needed
if ~isempty(subject.channels)
  cfg=[];
  cfg.channel = subject.channels;
  dataNoArtfct = ft_selectdata(cfg, dataNoArtfct);
end


%% Do ICA
if doica
  % resample
  cfg=[];
  cfg.resamplefs = 150;
  cfg.detrend = 'no';
  dataResample = ft_resampledata(cfg, dataNoArtfct);
  
  % ICA
  cfg=[];
  cfg.method          = 'fastica';
  cfg.channel         = 'MEG';
  cfg.fastica.numOfIC = 50;
  cfg.randomseed      = subjects.randomseed;
  comp = ft_componentanalysis(cfg, dataResample);
  
  % Inspect ICA components
  cfg = [];
  cfg.channel = {comp.label{1:10}}; % components to be plotted
  cfg.layout = 'CTF275_helmet.mat'; % specify the layout file that should be used for plotting
  cfg.compscale = 'local';
  ft_databrowser(cfg, comp);
  
  input('Please enter the ICA components in datainfo_erfosc. press any key to continue');
  subs_comp = input('Please also enter the to be rejected components here:');
else
  subs_comp = [subject.ecgcomp, subject.eyecomp];
end

%% Reject ICA components and artifacts
% decompose the original data as it was prior to downsampling to 150Hz
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_orig     = ft_componentanalysis(cfg, dataNoArtfct);
comp_orig.topo = comp.topo; % topographies should be the same for full rank (comp_orig) and rank deficient (comp, because only 50 components and 275 channels)

cfg=[];
cfg.component = subs_comp;
dataNoIca = ft_rejectcomponent(cfg, comp_orig, dataNoArtfct);

%% remove saccade artifacts
if ~existArtifact
  artifact_EOG_saccade = [];
  if ~exist('visAngleX') || ~exist('visAngleY')
    visAngleX = artfctdef.eyepos.visAngleX;
    visAngleY = artfctdef.eyepos.visAngleY;
  end
  % trials are marked as artifact if fixation is broken with more than
  % amount of visual degrees specified in visDegOffFixation (function
  % input argument)
  if ~exist('visDegFixDot')
    visDegFixDot =  artfctdef.eyepos.visDegFixDot;
  end
  % visDegOffFixation = visDegOffFixation + visDegFixDot; % relative to border of fixation dot
  % mark event as artifact if fixation from fixation point is broken with
  % more than x visual degree.
  artifact_EOG_saccade = break_fixation(data, visDegOffFixation);
  
  artfctdef.eyesaccade.artifact = artifact_EOG_saccade;
end

% reject eye movements from data
cfg = [];
cfg.artfctdef.eyesaccade.artifact = artfctdef.eyesaccade.artifact;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-1 3.75];
dataClean = ft_rejectartifact(cfg, dataNoIca);

cfg=[];
cfg.channel = 'HLC*';
headmotion = ft_selectdata(cfg, dataClean);

cfg=[];
cfg.channel = 'MEG';
dataClean = ft_selectdata(cfg, dataClean);

nTrials.postClean = length(dataClean.trial);
idxM = find(dataClean.trialinfo(:,5)>0 & dataClean.trialinfo(:,6)>0 & dataClean.trialinfo(:,2)==0);
nTrials.valid = length(idxM);


