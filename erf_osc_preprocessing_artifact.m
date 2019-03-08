function erf_osc_preprocessing_artifact(subj, isPilot, existArtifact, visDegOffFixation)
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

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(existArtifact)
    existArtifact = false;
end
if nargin<4 || isempty(visDegOffFixation)
    visDegOffFixation = 5;
end

% initiate diary
ft_diary('on')

%% Load data and define trials
erf_osc_datainfo; % load subject specific info.

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
cfg.trialfun = 'erf_osc_mytrialfun';
cfg.trialdef.prestim = min(cfg.logfile.log.realBaselineDuration, cfg.logfile.log.setBaselineDuration);
cfg.trialdef.poststim = cfg.logfile.log.completeDurationGrating;
cfg.catchtrial = cfg.logfile.log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);
cfg2=cfg; % use this for finding EOG artifacts.
nTrialsPreClean = size(cfg.trl,1);

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
    trlind=[];
    for i=1:size(dataRejVis.cfg.artfctdef.summary.artifact,1)
        trlind(i)=find(data.sampleinfo(:,1)==dataRejVis.cfg.artfctdef.summary.artifact(i,1));
    end;
    trlind = [trlind, subjects(subj).badtrials];
    if trlind
        disp(trlind);
        data.badtrials = trlind;
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
    [~, artifact_jump] = ft_artifact_zvalue(cfg, data);
    
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
    [~, artifact_muscle] = ft_artifact_zvalue(cfg, data);
    
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
    [~, artifact_EOG_blink] = ft_artifact_zvalue(cfg2);
   
    %% Find EOG saccade artifacts
    % this algorithm marks events as artifacts if eye movements are made
    % with more than 1 visual degree off the central fixation point.
    
    
    % pixel coords for 1920x1080 screen, upper left is (0,0)
    % these dimensions can be found in PHYSICAL.INI or your presentation/PTB
    % settings.
    screenRight=1919;
    screenLeft=0;
    screenTop=0;
    screenBottom=1079;
    
    % range of voltages and recording range as defined in FINAL.INI
    minVoltage=-5;
    maxVoltage=5;
    minRange=0;
    maxRange=1;
    
    totDist=780; % distance eye-screen in mm
    
    % if you want to calculate the visual degree per pixel, physical dimensions
    % of the pixels have to be known (so projected screen size).
    screenWidth_y=270;
    screenWidth_x=480; % has to be 'projected width' in mm (not actual width)
    
    % convert ADC-voltref back to GAZE position
    xGaze=[];
    yGaze=[];
    for iTrial=1:size(data.trial,2);% these are trials with right-left voluntary saccades
        
        
        idxChanHor=find(strcmp(data.label,'UADC005'));
        idxChanVer=find(strcmp(data.label,'UADC006'));
        
        voltageHor=data.trial{iTrial}(idxChanHor,:);
        voltageVer=data.trial{iTrial}(idxChanVer,:);
        
        % see eyelink1000_UserManual for formulas:
        RangeHor = (voltageHor-minVoltage)./(maxVoltage-minVoltage); % voltage range proportion
        S_h = RangeHor.*(maxRange-minRange)+minRange; % proportion of screen width or height
        
        R_v = (voltageVer-minVoltage)./(maxVoltage-minVoltage);
        S_v = R_v.*(maxRange-minRange)+minRange;
        
        xGaze{iTrial} = S_h.*(screenRight-screenLeft+1)+screenLeft;
        yGaze{iTrial} = S_v.*(screenBottom-screenTop+1)+screenTop;
        
        
        % gaze position on center of screen (where fixation dot is)
        Xgaze_origin=xGaze;
        Ygaze_origin=yGaze;
        
        xGaze{iTrial} = xGaze{iTrial} - round(screenRight/2);%%presentation/PTB assumes 1 more pixel than there actually is.
        %presenation x-position = screen position + 1 pixel to the right:
        %presentation x-pos(0) = 1919/2
        
        yGaze{iTrial} = yGaze{iTrial} - round(screenBottom/2);
    end
    % Xgaze and Ygaze now are the eye positions on the screen (relative to
    % fixation point). We can now convert them to visual angles (relative to
    % fixation point).
    
    % if you want to know the visual angle per pixel, calculate this
    %
    % X
    screenResX=screenRight+1;
    visAngX_rad = 2 * atan(screenWidth_x/2/totDist);
    visAngX_deg = visAngX_rad * (180/pi);
    visAngPerPixX = visAngX_deg / screenResX;
    
    dotSize = 20; % fixation dot is 20 pixels wide (see concentric_grating_experiment)
    visDegFixDot = visAngPerPixX * dotSize;
    
    % Y
    screenResY=screenBottom+1;
    visAngY_rad = 2 * atan(screenWidth_y/2/totDist);
    visAngY_deg = visAngY_rad * (180/pi);
    visAngPerPixY = visAngY_deg / screenResY;
    %}
    
    % visual angle in reference to middle of the screen
    for iTrial = 1:size(data.trial,2)
        visAngleX{iTrial} = rad2deg(atan((xGaze{iTrial}/2)/(totDist*screenResX/screenWidth_x)))*2;
        visAngleY{iTrial} = rad2deg(atan((yGaze{iTrial}/2)/(totDist*screenResY/screenWidth_y)))*2;
    end
    
    % NOTE: the rest of the eye movement artifact detection is done after
    % ICA based on visAngleX and Y.
    
    %% Save artifacts
    artfctdef.badtrial = trlind;
    artfctdef.jump.artifact = artifact_jump;
    artfctdef.muscle.artifact = artifact_muscle;
    artfctdef.eyeblink.artifact = artifact_EOG_blink;
    artfctdef.eyepos.visAngleX = visAngleX;
    artfctdef.eyepos.visAngleY = visAngleY; 
    artfctdef.eyepos.visDegFixDot = visDegFixDot; % save eye position in artifact definition in order to find artifacts in it later.
    %     artfctdef.eyesaccade.artifact = artifact_EOG_saccade; % save
    %     saccade artifacts later
    if isPilot
        if ~exist(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/', subj), 'dir');
            mkdir(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/', subj));
        end
        save(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
    else
        if ~exist(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/', subj), 'dir');
            mkdir(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/', subj));
        end
        save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
    end
    
    
    %% Load artifacts (only if specified)
    % if artifacts were selected before, load them.
else
    if isPilot
        load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
    else
        load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
    end
    if ~exist('artfctdef.eyesaccade')
        sprintf('WARNING: possibly eye saccade artifacts have not yet been marked. These do not appear in artfctdef.')
        input('press a key to continue')
    end
end

%% Reject artifacts from complete data
cfg = [];
cfg.artfctdef = artfctdef;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-1 3.75];
dataNoArtfct = ft_rejectartifact(cfg, data);

% reject extra channels if needed
if ~isempty(subjects(subj).channels) 
    cfg=[];
    cfg.channel = subjects(subj).channels;
    dataNoArtfct = ft_selectdata(cfg, dataNoArtfct);
end
    
    

%% Compute ICA
if ~existArtifact
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
    comp = ft_componentanalysis(cfg, dataResample);
    
    % save
    if isPilot
        filename = sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_icaComp.mat', subj, subj);
    else
        filename = sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_icaComp.mat', subj, subj);
    end
    save(filename, 'comp')
    
    %% Inspect ICA components
    % load the ica components computed with erf_osc_preprocessing_ica
    
    cfg = [];
    cfg.channel = {comp.label{1:10}}; % components to be plotted
    cfg.layout = 'CTF275_helmet.mat'; % specify the layout file that should be used for plotting
    cfg.compscale = 'local';
    ft_databrowser(cfg, comp);
    
    input('Please enter the ICA components in datainfo_erf_osc. press any key to continue');
    
else
    if isPilot
        filename = sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_icaComp.mat', subj, subj);
    else
        filename = sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_icaComp.mat', subj, subj);
    end
    load(filename, 'comp')
end

%% Reject ICA components and artifacts
% decompose the original data as it was prior to downsampling to 150Hz
cfg           = [];
cfg.unmixing  = comp.unmixing;
cfg.topolabel = comp.topolabel;
comp_orig     = ft_componentanalysis(cfg, dataNoArtfct);
comp_orig.topo = comp.topo; % topographies should be the same for full rank (comp_orig) and rank deficient (comp, because only 50 components and 275 channels)

% JM: voor de zekerheid kijken hoe comp_orig.topo eruitziet
% vergelijken met comp.topo;
% als niet hetzelfde, comp_orig.topo = comp.topo;


erf_osc_datainfo; % load subject specific data for the latest version.
if isPilot
    subs_comp = [pilotsubjects(subj).ecgcomp, pilotsubjects(subj).eyecomp];
else
    subs_comp = [subjects(subj).ecgcomp, subjects(subj).eyecomp];
end
cfg=[];
cfg.component = subs_comp;
dataNoIca = ft_rejectcomponent(cfg, comp_orig, dataNoArtfct);

% save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/dataHalfClean.mat', subj), 'dataNoIca', '-v7.3')

%% remove saccade artifacts
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
for iTrial = 1:size(data.trial,2)
    isFixation = sqrt(visAngleX{iTrial}.^2 + visAngleY{iTrial}.^2) < visDegOffFixation; % note all the samples where fixation is good/bad
    %         samplenumber = data.sampleinfo(iTrial,1);
    %         while samplenumber <= data.sampleinfo(iTrial,2)
    samplenumber = 1;
    while samplenumber <= size(isFixation,2) % for the whole trial (in sample numbers)
        % if there is a new off-fixation event, note the first off-fixation sample
        if any(isFixation(1, samplenumber:end)==0)
            sampleStartOffFixation = (samplenumber -1) + find(isFixation(1, samplenumber:end)==0,1);
            artifact_EOG_saccade(end+1,1) = data.sampleinfo(iTrial,1) + sampleStartOffFixation;
            %                 samplenumber = sampleStartOffFixation+1;
            %if there is fixation returns after off-fixation, note last
            % off-fixation sample of that event
            if any(isFixation(1, sampleStartOffFixation + 1 : end)==1)
                sampleStopOffFixation = sampleStartOffFixation + find(isFixation(1, sampleStartOffFixation + 1 : end) ==1, 1) -1; % stop off-fixation is 1 sample before on fixation again.
                artifact_EOG_saccade(end,2) = data.sampleinfo(iTrial,1) + sampleStopOffFixation;
                samplenumber = sampleStopOffFixation +1;
break
            else
                artifact_EOG_saccade(end,2) = data.sampleinfo(iTrial,2);
                samplenumber = data.sampleinfo(iTrial,2);
            end
        else
            break
        end
    end
end

% save artifacts in existing artifact definition.
if isPilot
    load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
else
    load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
end
artfctdef.eyesaccade.artifact = artifact_EOG_saccade;
if isPilot
    save(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
else
    save(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_artifact.mat', subj, subj), 'artfctdef');
end

% reject eye movements from data
cfg = [];
cfg.artfctdef.eyesaccade.artifact = artifact_EOG_saccade;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-1 3.75];
dataClean = ft_rejectartifact(cfg, dataNoIca);

cfg=[];
cfg.channel = 'HLC*';
headmotion = ft_selectdata(cfg, dataClean);

cfg=[];
cfg.channel = 'MEG';
dataClean = ft_selectdata(cfg, dataClean);
% cfg=[];
% cfg.detrend = 'no';
% cfg.demean = 'no';
% cfg.resamplefs = 600;
% dataClean = ft_resampledata(cfg, dataClean);

nTrialsPostClean = length(dataClean.trial);
idxM = find(dataClean.trialinfo(:,5)>0 & dataClean.trialinfo(:,6)>0 & dataClean.trialinfo(:,2)==0);
nTrialsValid = length(idxM);

%% save
if isPilot
    filename = sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/sub-%03d_cleandata', subj, subj);
else
    filename = sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata', subj, subj);
end
save(fullfile([filename '.mat']), 'dataClean','nTrialsPreClean','nTrialsPostClean','nTrialsValid', '-v7.3');
save(fullfile([filename(1:end-9), 'headmotion', '.mat']), 'headmotion');
ft_diary('off')

end


