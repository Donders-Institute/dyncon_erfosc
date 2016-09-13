function erf_osc_preprocessing_artifact(subj, existArtifact, isPilot)
% This function is interactive and can not be automated and/or
% paralellized. Before calling this function, make sure you performed ICA
% decomposition on the data in order to select and reject ECG artifacts
% (erf_osc_preprocessing_ica).
% This function loads data from the ERF-oscillation experiment and segments
% it into trials. Next, by visual inspection bad channels/trials are
% rejected and jump, muscle and EOG artifacts are selected. ICA components
% belonging to cardiac activity and artifacts are rejected from the data.
% The data is then downsampled and saved.
%
% use:
% subj, specify subject number (default = 1)
% existArtifact, true if artficats were already selected and saved (default
% = true)
% isPilot, true if datasets belong to pilot subjects (default = true)


if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    existArtifact = false;
end
if isempty(existArtifact)
    existArtifact = false;
end
if nargin<3
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

sprintf('Make sure you performed ICA decomposition with erf_osc_preprocessing_ica')

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
    for i=1:length(dataRejVis.cfg.artfctdef.summary.artifact)
        trlind(i)=find(data.sampleinfo(:,1)==dataRejVis.cfg.artfctdef.summary.artifact(i));
    end;
    if trlind
        disp(trlind);
        data.badtrials = trlind;
    else
        data.badtrials = [];
    end
    %% Find jump artifacts
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
    
    [~, artifact_muscle] = ft_artifact_zvalue(cfg, data);
    
    %% Find EOG saccade artifacts
    % NOTE: this algorhythm marks artifacts if saccades are made with a
    % visual angle larger than 1 degree. It does not mark trials as
    % artifacts if eye gaze is off fixation per se.
    
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
    %{
% X
screenResX=screenRight+1;
visAngX_rad = 2 * atan(screenWidth_x/2/totDist);
visAngX_deg = visAngX_rad * (180/pi);
visAngPerPixX = visAngX_deg / screenResX;

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
    
    % find onset and offset latencies of sacades with ft_detect_movement
    cfg=[];
    cfg.channel = {'UADC005', 'UADC006'};
    dataEye = ft_selectdata(cfg, data);
    cfg                     = [];
    cfg.method              = 'velocity2D';
    cfg.velocity2D.mindur   = 14; % minimum microsaccade durantion in samples (default = 3) (=12ms Engbert & Kliegl)
    cfg.velocity2D.velthres = 6; % threshold for velocity outlier detection (default = 6);
    [~, movement]           = ft_detect_movement(cfg,dataEye);
    
    % movement contains saccade on- and offset samples. find the corresponding
    % visual angles and compute the difference (for x and y). save these and
    % select trial as containing eyemovement if visualangle>1
    trlIdx = zeros(1,size(movement, 1));
    saccadeAngle = zeros(2,size(movement, 1));
    artfctdef.eye.saccade.artifact = [];
    
    for iSaccade = 1:size(movement, 1)
        % find the trial where the saccade belongs to
        trlIdx(iSaccade) = find(dataEye.sampleinfo(:,1)<=movement(iSaccade,1) & dataEye.sampleinfo(:,2)>=movement(iSaccade,2));
        % how many samples after the start of the trial does the saccade take
        % place and when does it end?
        saccadeStart = movement(iSaccade,1) - dataEye.sampleinfo(trlIdx(iSaccade),1) + 1;
        saccadeEnd = movement(iSaccade,2) - dataEye.sampleinfo(trlIdx(iSaccade),1) + 1;
        % The difference in the visual angle w.r.t. fixation at the start and
        % end of the saccade is the saccade visual angle (column 1 horizontal,
        % column 2 vertical).
        saccadeAngle(1, iSaccade) = abs(visAngleX{trlIdx(iSaccade)}(saccadeEnd) - visAngleX{trlIdx(iSaccade)}(saccadeStart));
        saccadeAngle(2, iSaccade) = abs(visAngleY{trlIdx(iSaccade)}(saccadeEnd) - visAngleY{trlIdx(iSaccade)}(saccadeStart));
        if saccadeAngle(:, iSaccade)>1 % otherwise it is a microsaccade.
            artifact_EOG_saccade(end+1,:) = movement(iSaccade, 1:2);
        end
    end
    
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
    cfg2.artfctdef.zvalue.bpfreq     = [1 15];
    cfg2.artfctdef.zvalue.bpfiltord  = 4;
    cfg2.artfctdef.zvalue.hilbert    = 'yes';
    % feedback
    cfg2.artfctdef.zvalue.interactive = 'yes';
    
    [~, artifact_EOG_blink] = ft_artifact_zvalue(cfg2);
    
    %% Save artifacts
    artfctdef.eyeblink.artifact = artifact_EOG_blink;
    artfctdef.eyesaccade.artifact = artifact_EOG_saccade;
    artfctdef.jump.artifact = artifact_jump;
    artfctdef.muscle.artifact = artifact_muscle;
    artfctdef.badtrial = trlind;
    if isPilot
        save(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/pilot/%02d/artifact.mat', subj), 'artfctdef');
    else
        save(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/experiment/%02d/artifact.mat', subj), 'artfctdef');
    end
    
    
    %% Load artifacts (only if specified)
    % if artifacts were selected before, load them.
else
    if isPilot
        load(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/pilot/%02d/artifact.mat', subj), 'artfctdef');
    else
        load(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/experiment/%02d/artifact.mat', subj), 'artfctdef');
    end
end

%% Inspect ICA components
%
% load the ica components computed with erf_osc_preprocessing_ica
if isPilot
    load(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/pilot/%02d/icaComp.mat', subj), 'comp');
else
    load(sprintf('/home/electromag/matves/Data/ERF_oscillation/artifacts/experiment/%02d/icaComp.mat', subj), 'comp');
end

if ~existArtifact
    cfg = [];
    cfg.channel = {comp.label{1:10}}; % components to be plotted
    cfg.layout = 'CTF275.lay'; % specify the layout file that should be used for plotting
    cfg.compscale = 'local';
    ft_databrowser(cfg, comp);
    
    input('Please enter the ICA components in datainfo_erf_osc. press any key to continue');
end

%% Reject ICA components and artifacts
erf_osc_datainfo; % load subject specific data for the latest version.
if isPilot
    subs_comp = pilotsubjects(subj).icacomp;
else
    subs_comp = subjects(subj).icacomp;
end
cfg=[];
cfg.component = subs_comp;
dataNoEcg = ft_rejectcomponent(cfg, comp, data);
% end
%}
%% Reject artifacts from complete data
cfg = [];
cfg.artfctdef = artfctdef;
cfg.artfctdef.reject = 'complete'; % remove complete trials
cfg.artfctdef.crittoilim = [-0.5 3.75];
dataClean = ft_rejectartifact(cfg, dataNoEcg);

% save
if isPilot
    save(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean', '-v7.3');
else
    save(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean', '-v7.3');
end

end


