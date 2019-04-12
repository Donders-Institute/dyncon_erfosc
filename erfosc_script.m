% the rationale for the stuff done in this script is as follows:
% subject specific data handling is done outside functions
% algorithmic pipeline steps are performed inside functions
if ~exist('dogroupanalysis', 'var'), dogroupanalysis = false; end
if ~dogroupanalysis
  if ~exist('subj', 'var')
    error('a subject needs to be specified');
  end
  if ~exist('isPilot', 'var')
    fprintf('assuming that a regular subject is requested\n');
    isPilot = false;
  end
  if isnumeric(subj) && ~isPilot
    erfosc_datainfo;
    subject = subjects(subj);
  elseif isnumeric(subj) && isPilot
    subject = pilotsubjects(subj);
  elseif ~isnumeric(subj)
    error('a subject should be identified by means of a scalar number, otherwise the data handling won''t work');
  end
end
if ~exist('dodics_gamma', 'var'), dodics = false; end
if ~exist('dofreq', 'var'),       dofreq = false;       end
if ~exist('docomp', 'var'),       docomp = false;       end
if ~exist('getdata', 'var'),      getdata = false;      end
if ~exist('dofreq_short', 'var'), dofreq_short = false; end
if ~exist('dolcmv_norm', 'var'),  dolcmv_norm = false; end
if ~exist('docorrpow_lcmv', 'var'), docorrpow_lcmv = false; end
if ~exist('docorr_gamma_rt', 'var'), docorr_gamma_rt = false; end
if ~exist('dostat_pow_erf', 'var'), dostat_pow_erf = false; end
if ~exist('dopartialcorr', 'var'), dopartialcorr = false; end
if ~exist('dolowfreq', 'var'), dolowfreq = false; end
if ~exist('dohighfreq', 'var'), dohighfreq = false; end
if ~exist('dotfa', 'var'), dotfa = false; end

if docorrpow_lcmv, dodics = true; dolcmv_parc = true; dofreq_short = true; end
if docorr_gamma_rt, dofreq_short = true; end
if dodics,       dofreq = true; end
if dofreq,       getdata = true; end
if docomp,       getdata = true; end
if dofreq_short, getdata = true; end
if dolcmv_parc,  getdata = true; end
if dolcmv_parc,  dolcmv_norm = true; end
if ~dolowfreq && ~dohighfreq, dohighfreq = 1; end


% load some things for source analysis
load(fullfile(subject.mridir,'preproc','headmodel.mat'));
cfg=[];
cfg.comment = 'load subject`s headmodel';
headmodel = ft_annotate(cfg, headmodel);
load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
cfg.comment = 'load subject`s 2D headmodel';
sourcemodel = ft_annotate(cfg, sourcemodel);
load('atlas_subparc374_8k.mat');
cfg.comment = 'load 2D atlas with 374 parcels, based on the Conte 69 atlas';
atlas = ft_annotate(cfg, atlas);
load('cortex_inflated_shifted');
cfg.comment = 'load cortex template';
ctx = ft_annotate(cfg, ctx);

if ~dogroupanalysis
  % this chunk creates 2 data structures [-0.75 0.5]
  if getdata
    basedir       = [project_dir, sprintf('processed/sub-%03d/ses-meg01', subj))];
    filename_data = fullfile(basedir, sprintf('sub-%03d_cleandata.mat', subj));
    load(filename_data);
    
    cfg=[];
    cfg.comment = 'load preprocessed MEG data';
    dataClean = ft_annotate(cfg, dataClean);
    [data_onset, data_shift, data_resp] = erfosc_getdata(dataClean);
    clear dataClean;
  end
  
  if dotfa
    cfg                 = [];
    cfg.method          = 'template';
    cfg.template        = 'CTF275_neighb.mat';
    cfg.neighbours      = ft_prepare_neighbours(cfg, data_onset);
    cfg.method          = 'sincos';
    data_planar         = ft_megplanar(cfg, data_onset);
    
    cfg                 = [];
    cfg.output          = 'pow';
    cfg.channel         = 'MEG';
    cfg.method          = 'mtmconvol';
    cfg.toi             = -1.75:0.05:1.75;
    cfg.keeptrials      = 'yes';
    cfg.pad             = 6;
    
    if dolowfreq
      cfg.taper        = 'hanning';
      cfg.foi          = 2:2:30;% analysis 2 to 30 Hz in steps of 2 Hz
      cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    elseif dohighfreq
      cfg.taper     = 'dpss';
      cfg.tapsmofrq = 8; % 8 Hz freq smoothing on both sides
      cfg.foi       = 28:4:100;
      cfg.t_ftimwin = ones(length(cfg.foi),1).*(1/4);
    end
    tfa_planar   = ft_freqanalysis(cfg, data_planar);
    cfg=[];
    cfg.method = 'sum';
    tfa          = ft_combineplanar(cfg, tfa_planar_low);
    cfg=[];
    cfg.latency = [-1 -0.25];
    tfa_baseline = ft_selectdata(cfg, tfa_high);
    if dolowfreq
        tmpstr = 'low';
    elseif dohighfreq
        tmpstr = 'high';
    end
    datadir = [results_dir sprintf('%03d', subj)];
    save([datadir, sprintf('tfa_%s_onset', subj, tmpstr)], 'tfa', 'tfa_baseline')
    clear tfa_planar data_planar
  end
  
  % this chunk does spectral decomposition - needed for computing spatial
  % filters
  if dofreq
    if dolowfreq
      if ~exist('peakfreq', 'var')
        load([project_dir, sprintf('analysis/freq/sub-%03d/sub-%03d_pow_low.mat',subj,subj)]);
      end
      foi = [peakfreq peakfreq]; smo = 2;
      latency = [-0.75+1./data_onset.fsample 0-1./data_onset.fsample];
      [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    elseif dohighfreq
      if ~exist('latency', 'var')
        latency = [-inf 0-1./data_onset.fsample];
      end
      [freq_onset, freq_shift, P, pow_onset, pow_shift] = erfosc_freq(data_onset, data_shift, latency, subject);
    end
    if dosave
      savedir = [project_dir 'analysis/freq'];
      save(fullfile(savedir, sprintf('sub-%03d/sub-%03d_freq', subj,subj)), 'freq_shift', 'P', 'latency');
    end
  end
  
  % this chunk does source reconstruction
  if dodics
    [sourcefreq_onset, sourcefreq_shift, sourcefreq_shift_Tval] = erfosc_dics(freq_onset, freq_shift, headmodel, sourcemodel);
      if dolowfreq
        tmpstr = 'low';
      elseif dohighfreq
        tmpstr = 'high';
      end
      savedir = [results_dir, sprintf('%03d', subj)];
      filename = fullfile(savedir, sprintf('%spow_tval.mat', tmpstr));
      save(filename, 'sourcefreq_shift_Tval');
  end
  
  % this chunk does spectral decomposition - used for estimating
  % gamma/lowfreq power before ERF peak
  if dofreq_short
    if dohighfreq
      latency = subject.erflatency(1).*[1 1] - 0.02 - [0.20 1./600];%MvE
      foi     = subject.gammapeak(end).*[1 1];
      smo     = max(10, diff(subject.gammaband(end,:))./2);
    elseif dolowfreq
      load([project_dir, sprintf('analysis/freq/sub-%03d/sub-%03d_pow_low.mat',subj,subj)]);
      latency = subject.erflatency(1).*[1 1] - 0.02 - [0.4 1./600];%MvE
      foi = peakFreq;
      smo     = 2.5;
    end
    foi = foi.*[1 1];
    [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    if ~exist('savefreq', 'var')
      savefreq = false;
    end
    if savefreq
      if dolowfreq
        tmpstr = 'low';
      elseif dohighfreq
        tmpstr = 'high';
      end
      savedir = [project_dir 'analysis/freq/'];
      filename = fullfile(savedir, sprintf('sub-%03d/sub-%03d_%sfreqshort', subj,subj, tmpstr));
      save(filename, 'freq_shift', 'P', 'latency');
    end
  end
  
  if dodics && dofreq_short
    sourcefreq_shift.pow = transpose((abs(sourcefreq_shift.F*transpose(freq_shift.fourierspctrm)).^2)*P);
    cfg=[];
    cfg.comment = 'for pow take the power in the pre-ERF window by combining the spatial filter with the fourierspectrum of this window';
    sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);
    sourcefreq_shift.pow = standardise(log10(sourcefreq_shift.pow), 1);
    cfg.comment = 'log10 avg.pow and standardise w.r.t. trials';
    sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);
  end
  
  if dolcmv_parc
    [source_parc] = erfosc_lcmv_parc(data_resp, headmodel, sourcemodel, atlas);
    if dosave
      savedir = [project_dir 'analysis/source/'];
      save(fullfile(savedir, sprintf('sub-%03d/sub-%03d_lcmv',subj, subj)), 'source_parc');
    end
  end
  
  if dolcmv_norm
    tmpcfg = [];
    tmpcfg.latency = [-0.2 -1./600];
    data_shift_baseline = ft_selectdata(tmpcfg, data_shift);
    
    tmpcfg = [];
    tmpcfg.covariance = 'yes';
    tmp = ft_timelockanalysis(tmpcfg, data_shift_baseline);
    
    source_parc.noise = diag(sqrt(source_parc.F*tmp.cov*source_parc.F'));
    cfg=[];
    cfg.comment = 'compute the noise by calculating "noise = diag(sqrt(lcmv_spatial_filters*noise_covariance*transpose(lcmv_spatial_filters)))"'; %FIXME
    source_parc = ft_annotate(cfg, source_parc);
    
    if doplot && subj==13
      cfg = [];
      cfg.preproc.demean='yes';
      cfg.preproc.baselinewindow = [-.1 0];
      cfg.vartrllength = 2;
      tmptlck = ft_timelockanalysis(cfg, data_shift);
      cfg=[];
      cfg.channel = '';
    end
    
    if dosave
      save(filename, 'noise', '-append');
    end
  end
  
  if docorrpow_lcmv
    cfg=[];
    cfg.matrix = repmat(1./source_parc.noise, [1 size(source_parc.avg,2)]);
    cfg.parameter = 'avg';
    cfg.operation = 'multiply';
    cfg.comment = 'devide the source level timelock average by the noise';
    source_parc_norm = ft_math(cfg, source_parc);
    source_parc_norm = copyfields(source_parc, source_parc_norm, {'noise', 'F'});
    cfg=[];
    cfg.comment = 'copy fields noise and spatial filter (F) from before ft_math';
    source_parc_norm = ft_annotate(cfg, source_parc_norm);
    
    cfg=[];
    cfg.baseline = [-0.1 -1./data_shift.fsample];
    source_parc_norm = ft_timelockbaseline(cfg, source_parc_norm);
    
    cfg = [];
    cfg.latency = [-0.1 0.5-1./data_shift.fsample];
    datapst = ft_selectdata(cfg, data_shift);
    
    datapst.trial = source_parc_norm.F*datapst.trial;
    datapst.label = source_parc_norm.label;
    cfg=[];
    cfg.comment = 'mutliply the LCMV spatial filter with the single-trial timelock data.';
    datapst = ft_annotate(cfg, datapst);
    
    cfg = [];
    cfg.demean = 'yes';
    cfg.baselinewindow = [-inf 0];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    cfg.lpfilttype = 'firws';
    cfg.lpfiltdir = 'onepass-reverse-zerophase'
    datapst = ft_preprocessing(cfg, datapst);
    
    cfg=[];
    cfg.keeptrials = 'yes';
    datapst = ft_timelockanalysis(cfg, datapst);
    cfg = [];
    cfg.latency = subject.erflatency;
    cfg.avgovertime = 'yes';
    datapeak = ft_selectdata(cfg,datapst);
    
    X = datapeak.trial;
    Xpow = abs(mean(X,1));
    signswap = repmat(sign(mean(X,1)), [size(datapeak.trial,1), 1]);
    
    % let the amplitude be on average positive
    % FIXME when issue 1068 has been resolved, replace this:
    datapeak.trial = datapeak.trial.*signswap;
    cfg=[];
    cfg.comment = 'multiply the peak data by a matrix which lets the trial-average be positive at every parcel.';
    datapeak = ft_annotate(cfg, datapeak);
    % by this:
    %{
  
  cfg=[];
  cfg.operation = 'multiply';
  cfg.matrix = signswap;
  cfg.parameter = 'trial';
  datapeak = ft_math(cfg, datapeak);
    %}
    datapeak.trial = standardise(datapeak.trial,1);
    cfg=[];
    cfg.comment = 'Standardise peak data over trials';
    datapeak = ft_annotate(cfg, datapeak);
    
    cfg = [];
    datapst = ft_timelockanalysis(cfg, datapst);
    cfg = [];
    cfg.matrix = repmat(signswap(1,:), [size(datapst.avg,2), 1])';
    cfg.operation = 'multiply';
    cfg.comment = sprintf('align the parcels such that the trial-average of the [%d %d]s window is positive.', subject.erflatency);
    cfg.parameter = 'avg';
    datapst = ft_math(cfg, datapst);
    % do the correlations
    
    if dolowfreq
      [m, idx] = max(sourcefreq_shift_Tval.stat);
      cfg=[];
      cfg.comment = 'Create field roipow: from avg.pow, select only single-trial power values at the source location of maximal NEGATIVE induced power change';
    elseif dohighfreq
      [m, idx] = max(sourcefreq_shift_Tval.stat);
      cfg=[];
      cfg.comment = 'Create field roipow: from avg.pow, select only single-trial power values at the source location of maximal POSITIVE induced power change';
    end
    sourcefreq_shift.roipow = sourcefreq_shift.pow(:,idx);
    sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);
    
    sourcepow = rmfield(sourcefreq_shift, {'F', 'avg'});
    sourcepow.trialinfo = datapeak.trialinfo;
    cfg=[];
    cfg.comment = 'create source structure with single-trial power on source level.';
    sourcepow = ft_annotate(cfg, sourcepow);
    
    sourcepow.rho_pow_rt = corr(sourcepow.pow, sourcepow.trialinfo(:,7), 'type', 'spearman');
    sourcepow.rho_roipow_rt = corr(sourcepow.roipow, sourcepow.trialinfo(:,7), 'type', 'spearman');
%     jitter = log10(sourcepow.trialinfo(:,5)-sourcepow.trialinfo(:,4));
%     jitter = (jitter-mean(jitter))./1200; % demean and divide by sampling frequency
%     sourcepow.rho_roipow_rt = corr(sourcepow.roipow, sourcepow.trialinfo(:,7), 'type', 'spearman');
%     sourcepow.partialrho_roipow_jitter = corr(sourcepow.roipow,  jitter, 'type', 'spearman');
    
    cfg=[]
    cfg.comment = 'calculate Spearman rank correlation between power and reaction times (7th column trialinfo).';
    datapeak = ft_annotate(cfg, datapeak);
    
    datapeak.rho_pow_erf = corr(datapeak.trial, sourcefreq_shift.roipow, 'type', 'spearman');
    cfg=[];
    cfg.comment = 'calculate Spearman rank correlation between roipow and peak ERF amplitude.';
    datapeak = ft_annotate(cfg, datapeak);
    
    if dopartialcorr % FIXME errors in ft_timelockanalysis and ft_regressconfound in new fieldtrip version.
      % prepare eye data
      [~, eyedata_shift] = erfosc_preprocessing_eyedata(subj, data_onset);
      % get gamma-ERF correlation, accounting for pupil diameter, without confounding eye position.
      cfg=[];
      cfg.vartrllength = 2;
      cfg.keeptrials = 'yes';
      eyedata_shift = ft_timelockanalysis(cfg,eyedata_shift);
      cfg=[];
      cfg.comment = 'get pupil diameter data without contribution of eye position';
      cfg.toilim = [-0.2 -1./600] + subject.erflatency(1) - 0.02;
      pupil_diameter = erfosc_regress_eye(ft_redefinetrial(cfg, eyedata_shift), {'UADC007'}, {'visAngleX', 'visAngleY'});
      cfg=[];
      cfg.avgovertime = 'yes';
      pupil_diameter = ft_selectdata(cfg, pupil_diameter);
      pupil_diameter.trial = standardise(pupil_diameter.trial);
      cfg=[];
      cfg.comment = 'standardise pupil diameter w.r.t. trials';
      pupild = ft_annotate(cfg, pupild);
      
      idx = match_str(eyedata_shift.label, {'visAngleX', 'visAngleY'});
      eyedata_shift.trial(:,end+1,:) = (eyedata_shift.trial(:,idx(1),:).^2 + eyedata_shift.trial(:,idx(2),:).^2).^0.5;
      eyedata_shift.label{end+1} = 'distance';
      cfg=[];
      cfg.comment = 'add channel `distance` as distance from fixation (combining vertical and horizontal distance).';
      eyedata_shift = ft_annotate(cfg, eyedata_shift);
      cfg=[];
      cfg.comment = 'get eye position data without contribution of pupil diameter';
      cfg.toilim = [-0.2 -1./600] + subject.erflatency(1) - 0.02;
      eyepos = erfosc_regress_eye(ft_redefinetrial(cfg, eyedata_shift), {'distance'}, {'UADC007'});
      cfg=[];
      cfg.avgovertime = 'yes';
      distance_fixation = ft_selectdata(cfg, eyepos);
      distance_fixation.trial = standardise(distance_fixation.trial);
      cfg=[];
      cfg.comment = 'standardise gaze distance from fixation w.r.t. trials';
      distance_fixation = ft_annotate(cfg, distance_fixation);
      
      % partial correlations with eye data
      datapeak.partialrho_pupild = partialcorr(datapeak.trial, sourcefreq_shift.roipow, pupil_diameter.trial, 'type', 'spearman');
      cfg=[];
      cfg.comment = 'add field `partialrho_pupild`, the partial Spearman correlation betweek roipow and ERF peak amplitude, accounted for pupil diameter.';
      datapeak = ft_annotate(cfg, datapeak);
      
      % get correlation gamma-ERF, accounting for eye position, without confound pupil diameter
      datapeak.partialrho_distancefixation = partialcorr(datapeak.trial, sourcefreq_shift.roipow, distance_fixation.trial, 'type', 'spearman');
      cfg=[];
      cfg.comment = 'add field `partialrho_distancefixation`, the partial Spearman correlation betweek roipow and ERF peak amplitude, accounted for gaze distance to fixation.';
      datapeak = ft_annotate(cfg, datapeak);
      % get correlation gamma-ERF, accounting for eye position and pupil
      % diameter (both not confounded by the other)
      datapeak.partialrho_pupuld_distancefixation = partialcorr(datapeak.trial, sourcefreq_shift.roipow, [distance_fixation.trial, pupil_diameter.trial], 'type', 'spearman');
      cfg=[];
      cfg.comment = 'add field `partialrho_pupild_distancefixation`, the partial Spearman correlation betweek roipow and ERF peak amplitude, accounted for pupil diameter and gaze distance to fixation (after regressing out their interaction).';
      datapeak = ft_annotate(cfg, datapeak);
    end
        
    datapeak.rho_erf_rt = corr(datapeak.trial, datapeak.trialinfo(:,7), 'type', 'spearman'); %MvE
    cfg=[]
    cfg.comment = 'calculate Spearman rank correlation between ERF amplitude and reaction times (7th column trialinfo).';
    datapeak = ft_annotate(cfg, datapeak);
    
    exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
    
    source = [];
    source.freq          = 0;
    source.brainordinate = atlas;
    source.label         = atlas.parcellationlabel;
    source.rho_pow_erf           = nan(374,1);
    source.rho_erf_rt            = nan(374,1);
    if dopartialcorr
      source.partialrho_pupild                  = nan(374,1);
      source.partialrho_distancefixation        = nan(374,1);
      source.partialrho_pupild_distancefixation = nan(374,1);
    end
    source.pow           = nan(374,1);
    source.dimord        = 'chan_freq';
    
    indx = 1:374;
    indx(exclude_label) = [];
    source.rho_pow_erf(indx)    = datapeak.rho_pow_erf;
    source.rho_erf_rt(indx)     = datapeak.rho_erf_rt;
    cfg=[];
    if dolowfreq
      tmpstr = 'low';
    elseif dohighfreq
      tmpstr = 'high';
    end
    if dopartialcorr
      source.partialrho_pupild = datapeak.partialrho_pupild;
      source.partialrho_distancefixation = datapeak.partialrho_distancefixation;
      source.partialrho_pupild_distancefixation = datapeak.partialrho_pupild_distancefixation;
      cfg.comment = sprintf('mold all correlations into a fieldtrip source structure, as well as source-level ERF and %s power. Contains rho_pow_erf, partialrho_pupild, partialrho_distancefixation, partialrho_pupild_distancefixation', tmpstr);
    else
      cfg.comment = sprintf('mold correlation into a fieldtrip source structure, as well as source-level ERF and %s power', tmpstr);
    end

    source.pow(indx)    = Xpow(:);
    source.trial        = datapeak.trial;
    source = ft_annotate(cfg, source);
    
    % save the results to disk.
    datadir = [results_dir sprintf('%03d', subj)];
    save(fullfile(datadir, sprintf('source_erf_%spow', tmpstr)), 'source');
    save(fullfile(datadir, sprintf('source_%spow', tmpstr)), 'sourcepow');
  end
else
  if dotfa
    k=1;
    for subj = allsubs
      datadir = [results_dir sprintf('%03d', subj)];
      if dolowfreq
        tmpstr = 'low';
      elseif dohighfreq
        tmpstr = 'high';
      end
      filename1 = fullfile(datadir, sprintf('tfa_%s_onset', tmpstr));
      filename2 = fullfile(datadir, sprintf('%spow_tval', tmpstr));
      fprintf('loading data for subject %s\n', subj);
      load(filename1, filename2);
      freq{k} = tfa;
      baseline{k} = tfa_baseline;
      pow_tval{k} = sourcefreq_shift_Tval;
      k=k+1;
    end
   
    cfg=[];
    cfg.comment = 'Compute relative change of TFR wrt mean baseline over time';
    for k=1:numel(freq)  
      baseline{k}.powspctrm = repmat(mean(baseline{k}.powspctrm,3),[1 1 numel(freq{k}.time)]);
      freq{k}.powspctrm = freq{k}.powspctrm./baseline{k}.powspctrm-1;
      freq{k} = ft_annote(cfg, freq{k});
    end
    
    cfg=[];
    cfg.appenddim = 'rpt';
    freq_allsubs = ft_appendfreq(cfg, freq{:});
    pow_tval_allsubs = ft_appendfreq(cfg, pow_tval{:});
    if dolowfreq
      [~, idx] = min(pow_tval_allsubs.stat,[], 2); % maximum tvalues
    elseif dohighfreq
      [~, idx] = max(pow_tval_allsubs.stat,[], 2); % maximum tvalues
    end
    cfg=[];
    cfg.avgoverrpt = 'yes';
    pow_tval_allsubs = ft_selectdata(cfg, pow_tval_allsubs);
    
    
    if doplot
      virtualchanpos=[];
      virtualchanpos.chanpos = ctx.pos(idx,:);
      virtualchanpos.elecpos = ctx.pos(idx,:);
      for k=1:numel(allsubs)
        virtualchanpos.label{k}   = sprintf('%d', k);
      end
      pow_tval_allsubs.pos = ctx.pos;
      pow_tval_allsubs.virtualchanpos = virtualchanpos;
      cfg=[];
      cfg.comment = 'replace the .pos info by the shifted version. add virtualchanpos info.';
      pow_tval_allsubs = ft_annotate(cfg, pow_tval_allsubs);
      plotpart = 1; erfosc_plotting;
      
      
    end
  end
  
  if dostat_pow_erf
    k=1;
    for subj = allsubs
      datadir = [results_dir sprintf('%03d', subj)];
      if dolowfreq
        tmpstr = 'low';
      elseif dohighfreq
        tmpstr = 'high';
      end
      filename = fullfile(datadir, sprintf('source_erf_%spow', tmpstr));
      fprintf('loading data for subject %s\n', subj);
      load(filename);
      S{k} = source;
      k=k+1;
    end
    
    n = 32;
    
    S0{1} = S{1};
    S0{1}.rho_pow_erf = 0*S0{1}.rho_pow_erf;
    S0{1}.rho_erf_rt  = 0*S0{1}.rho_erf_rt;
    if dopartialcorr
      S0{1}.partialrho_pupild = 0*S0{1}.partialrho_pupild;
      S0{1}.partialrho_distancefixation = 0*S0{1}.partialrho_distancefixation;
      S0{1}.partialrho_pupild_distancefixation = 0*S0{1}.partialrho_pupild_distancefixation;
    end
    cfg=[];
    cfg.comment = 'multiply rho with zero to create a zero-distribution.';
    S0{1} = ft_annotate(cfg, S0{1});
    for k=2:n
      S0{k} = S0{1};
    end
    
    tmp = full(parcellation2connectivity_midline(atlas));
    for k=1:size(tmp,1)
      neighb(k).label = atlas.parcellationlabel{k};
      neighb(k).neighblabel = atlas.parcellationlabel(find(tmp(k,:)));
    end
    
    cfgs = [];
    cfgs.comment = 'compare the correlations to zero';
    cfgs.method='montecarlo';
    cfgs.design=[ones(1,n) ones(1,n)*2;1:n 1:n];
    cfgs.statistic='ft_statfun_wilcoxon';
    cfgs.numrandomization=10000;
    cfgs.ivar=1;
    cfgs.uvar=2;
    cfgs.correctm='cluster';
    cfgs.clusterthreshold='nonparametric_individual';
    cfgs.connectivity = neighb;
    cfgs.neighbours = neighb;
    cfgs.tail = 1;
    cfgs.clustertail = 1;
    cfgs.correcttail = 'prob';
    cfgs.clusteralpha = 0.05;
    cfgs.alpha = 0.05;
    
    cfgs.parameter = 'rho_pow_erf';
    stat_pow_erf = ft_freqstatistics(cfgs,S{:},S0{:});
    
    cfgs.parameter = 'rho_erf_rt';
    stat_erf_rt = ft_freqstatistics(cfgs,S{:},S0{:});

    
    if dopartialcorr
      cfgs.parameter = 'partialrho_pupild';
      stat_pupild=ft_freqstatistics(cfgs,S{:},S0{:});
      cfgs.parameter = 'partialrho_distancefixation';
      stat_xy=ft_freqstatistics(cfgs,S{:},S0{:});
      cfgs.parameter = 'partialrho_pupild_distancefixation';
      stat_eye=ft_freqstatistics(cfgs,S{:},S0{:});
    end
    datadir = [results_dir 'group/'];
    filename = [project_dir, sprintf('stat_erf_%spow.mat', tmpstr)];
    if dopartialcorr
      save(filename, 'S', 'stat_pow_erf', 'stat_erf_rt', 'stat_eye', 'stat_pupild', 'stat_xy');
    else
      save(filename, 'S', 'stat_pow_erf', 'stat_erf_rt');
    end
  end
end