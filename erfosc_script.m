% the rationale for the stuff done in this script is as follows:
% subject specific data handling is done outside functions
% algorithmic pipeline steps are performed inside functions

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

if ~exist('dodics_gamma', 'var'), dodics_gamma = false; end
if ~exist('dodics_lowfreq', 'var'), dodics_lowfreq = false; end
if ~exist('dofreq', 'var'),       dofreq = false;       end
if ~exist('docomp', 'var'),       docomp = false;       end
if ~exist('getdata', 'var'),      getdata = false;      end
if ~exist('dofreq_short', 'var'), dofreq_short = false; end
if ~exist('dofreq_short_lowfreq', 'var'), dofreq_short_lowfreq = false; end
if ~exist('dolcmv_norm', 'var'),  dolcmv_norm = false; end
if ~exist('docorrpow_lcmv', 'var'), docorrpow_lcmv = false; end
if ~exist('docorrpow_lcmv_lowfreq', 'var'), docorrpow_lcmv_lowfreq = false; end
if ~exist('docorr_gamma_rt'), docorr_gamma_rt = false; end
if ~exist('dostat_erf_rt'), dostat_erf_rt = false; end
if ~exist('dostat_pow_erf'), dostat_pow_erf = false; end

if docorr_gamma_rt, dofreq_short = true; end
if dodics_gamma, dofreq  = true; end
if dodics_lowfreq, dofreq  = true; end
if dofreq,       getdata = true; end
if docomp,       getdata = true; end
if dofreq_short, getdata = true; end
if dofreq_short_lowfreq, getdata = true; end
if dolcmv_parc,  getdata = true; end
if dolcmv_norm,  getdata = true; end



% this chunk creates 2 data structures [-0.75 0.5]
if getdata
  if ~exist('undocomp','var')
    undocomp = false;
  end
  [p,f,e]       = fileparts(subject.dataset);
  basedir       = strrep(p, 'raw', 'processed');
  filename_data = fullfile(basedir, sprintf('sub-%03d_cleandata.mat', subj));
  load(filename_data);
  
  cfg=[];
  cfg.comment = 'load preprocessed MEG data';
  dataClean = ft_annotate(cfg, dataClean);
  if undocomp
    filename_comp = fullfile(basedir, sprintf('sub-%03d_icaComp.mat', subj));
    load(filename_comp);
    sel = zeros(0,1);
    if ~isempty(subject.ecgcomp)
      sel = [sel subject.ecgcomp];
    end
    if ~isempty(subject.eyecomp)
      sel = [sel subject.eyecomp];
    end
    comp.topo     = comp.topo(:,sel);
    comp.unmixing = comp.unmixing(sel,:);
  else
    comp = [];
  end
  [data_onset, data_shift, data_resp] = erfosc_getdata(dataClean, comp);
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

  savedir = [project_dir, sprintf('results/freq/sub-%03d', subj)];
  if ~exist('dolowfreq', 'var'); dolowfreq = 1; end
  if ~exist('dohighfreq', 'var'); dohighfreq = 1; end
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
  if dosave
    if dolowfreq
      tmpstr = 'low';
    elseif dohighfreq
      tmpstr = 'high';
    end
    save([savedir, sprintf('sub-%03d_tfa_%s_onset', subj, tmpstr)], 'tfa', 'tfa_baseline')
  end
  clear tfa_planar_low tfa_planar_high data_planar
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
  load(fullfile(subject.mridir,'preproc','headmodel.mat'));
  cfg=[];
  cfg.comment = 'load subject`s headmodel';
  headmodel = ft_annotate(cfg, headmodel);
  load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
  cfg.comment = 'load subject`s 2D headmodel';
  sourcemodel = ft_annotate(cfg, sourcemodel);
  [sourcefreq_onset, sourcefreq_shift, sourcefreq_shift_Tval] = erfosc_dics(freq_onset, freq_shift, headmodel, sourcemodel);
  
  savedir = [project_dir 'analysis/source/'];
  if dosave
    filename = fullfile(savedir, sprintf('sub-%03d/sub-%03d_source', subj,subj));
    if dolowfreq
      filename = [filename, '_low'];
    end
    save(filename, 'sourcefreq_onset', 'sourcefreq_shift', 'sourcefreq_shift_Tval', 'F');
  end
end

% this chunk does spectral decomposition - used for estimating
% gamma/lowfreq power before ERF peak
if dofreq_short
  peakpicking;
  if dohighfreq
    latency = peaks(subj,1).*[1 1] - 0.02 - [0.20 1./600];%MvE
    foi     = subject.gammapeak(end).*[1 1];
    smo     = max(10, diff(subject.gammaband(end,:))./2);
  elseif dolowfreq
    load([project_dir, sprintf('analysis/freq/sub-%03d/sub-%03d_pow_low.mat',subj,subj)]);
    latency = peaks(subj,1).*[1 1] - 0.02 - [0.4 1./600];%MvE
    foi = peakFreq;
    smo     = 2.5;
  end
  foi = foi.*[1 1];
  [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
  if ~exist('savefreq', 'var')
    savefreq = false;
  end
  if savefreq
    savedir = [project_dir 'analysis/freq/'];
    filename = fullfile(savedir, sprintf('sub-%03d/sub-%03d_freqshort', subj,subj));
    if dolowfreq
      filename = [filename, '_low'];
    end
    save(filename, 'freq_shift', 'P', 'latency');
  end
end

if dolcmv_parc
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
  if dosave
    save(filename, 'noise', '-append');
  end
end

if docorrpow_lcmv
  peakpicking;
  cfg=[];
  cfg.matrix = repmat(1./source_parc.noise, [1 size(source_parc.avg,2)]);
  cfg.parameter = 'avg';
  cfg.operation = 'multiply';
  cfg.comment = 'devide the average by the noise';
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
  cfg.comment = 'mutliply the LCMV spatial filter with the data (time locked to stimulus change)';
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
  cfg.latency = peaks(subj,:);
  cfg.avgovertime = 'yes';
  datapeak = ft_selectdata(cfg,datapst);
  
  X = datapeak.trial;
  Xpow = abs(mean(X,1));
  signswap = repmat(sign(mean(X,1)), [size(datapeak.trial,1), 1]);
  
  % let the amplitude be on average positive
  %{
  % FIXME when issue 1068 has been resolved, replace this:
  datapeak.trial = datapeak.trial.*signswap;

  cfg=[];
  cfg.comment = 'multiply the peak data by a matrix which lets the trial-average be positive at every parcel.';
  datapeak = ft_annotate(cfg, datapeak);
  % by this:
  %}
  
  cfg=[];
  cfg.operation = 'multiply';
  cfg.matrix = signswap;
  cfg.parameter = 'trial';
  datapeak = ft_math(cfg, datapeak);
  
  datapeak.trial = standardise(datapeak.trial,1);
  cfg=[];
  cfg.comment = 'Standardise peak data over trials';
  datapeak = ft_annotate(cfg, datapeak);
  
  cfg = [];
  datapst = ft_timelockanalysis(cfg, datapst);
  cfg = [];
  cfg.matrix = repmat(signswap(1,:), [size(datapst.avg,2), 1])';
  cfg.operation = 'multiply';
  cfg.comment = sprintf('align the parcels such that the trial-average of the [%d %d]s window is positive.',peaks(subj,:));
  cfg.parameter = 'avg';
  datapst = ft_math(cfg, datapst);
  
  % do the correlations
  if dolowfreq
    [m, idx] = max(sourcefreq_shift_Tval.stat);
    sourcefreq_shift.roipow      = (abs(sourcefreq_shift.F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    cfg=[];
    cfg.comment = 'add `roipow` field to frequency structure, with single-trial power values at the source location of maximal NEGATIVE induced power change';
  elseif dohighfreq
    % prepare eye data
    [~, eyedata_shift] = erfosc_preprocessing_eyedata(subj);
    % get gamma-ERF correlation, accounting for pupil diameter, without confounding eye position.
    cfg=[];
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
    eyedata_shift = ft_timelockanalysis(cfg,eyedata_shift);
    cfg=[];
    cfg.comment = 'get pupil diameter data without contribution of eye position';
    cfg.toilim = [-0.2 -1./600] + peaks(subj,1) - 0.02;
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
    cfg.toilim = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    eyepos = erfosc_regress_eye(ft_redefinetrial(cfg, eyedata_shift), {'distance'}, {'UADC007'});
    cfg=[];
    cfg.avgovertime = 'yes';
    distance_fixation = ft_selectdata(cfg, eyepos);
    distance_fixation.trial = standardise(distance_fixation.trial);
    cfg=[];
    cfg.comment = 'standardise gaze distance from fixation w.r.t. trials';
    distance_fixation = ft_annotate(cfg, distance_fixation);
    
    [m, idx] = max(sourcefreq_shift_Tval.stat);
    sourcefreq_shift.roipow      = (abs(sourcefreq_shift.F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    cfg=[];
    cfg.comment = 'add `roipow` field to frequency structure, with single-trial power values at the source location of maximal POSITIVE induced power change';
  end
  sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);
  cfg=[];
  cfg.parameter = 'roipow';
  cfg.operation = 'log10';
  sourcefreq_shift = ft_math(cfg, sourcefreq_shift);
  
  sourcefreq_shift.roipow      = standardise(sourcefreq_shift.roipow(:));
  cfg=[];
  cfg.comment = 'standardise roipow w.r.t. trials';
  sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);
  
  datapeak.rho = corr(datapeak.trial, sourcefreq_shift.roipow, 'type', 'spearman');
  cfg=[];
  cfg.comment = 'calculate Spearman rank correlation between roipow and peak ERF amplitude.';
  datapeak = ft_annotate(cfg, datapeak);
  
  if dohighfreq
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
  
  load atlas_subparc374_8k
  cfg=[];
  cfg.comment = 'load 2D atlas with 374 parcels, based on the Conte 69 atlas';
  atlas = ft_annotate(cfg, atlas);
  exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
  
  source = [];
  source.brainordinate = atlas;
  source.label         = atlas.parcellationlabel;
  source.rho           = zeros(374,1);
  source.partialrho    = zeros(374,1);
  source.pow           = zeros(374,1);
  source.dimord        = 'chan';
  
  indx = 1:374;
  indx(exclude_label) = [];
  source.rho(indx)    = datapeak.rho;
  cfg=[];
  if dohighfreq
    source.partialrho(indx,1) = datapeak.partialrho_pupild;
    source.partialrho(indx,2) = datapeak.partialrho_distancefixation;
    source.partialrho(indx,3) = datapeak.partialrho_pupild_distancefixation;
    cfg.comment = 'mold all correlations into a fieldtrip source structure. The field `partialrho` contains partialrho_pupild, partialrho_distancefixation, partialrho_pupild_distancefixation in columns 1, 2, and 3, respectively.'
    tmpstr = 'low';
  else
    cfg.comment = 'mold all correlations into a fieldtrip source structure.';
    tmpstr = 'gamma';
  end
  source.pow(indx)    = Xpow(:);
  cfg=[];
  cfg.comment = 'mold all correlations into a fieldtrip source structure. The field `partialrho` contains partialrho_pupild, partialrho_distancefixation, partialrho_pupild_distancefixation in columns 1, 2, and 3, respectively.'
  source = ft_annotate(cfg, source);
  if dosave
    datadir = [project_dir 'analysis/corr/'];
    save(fullfile(datadir, sprintf('sub-%03d/sub-%03d_corrpowlcmv_%s',subj,  subj, tmpstr)), 'source', 'sourcefreq_shift', 'datapeak')%, 'pupild', 'distance');
  end
end

if docorr_gamma_rt
  datadir = [project_dir 'analysis/'];
  load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_source',  subj,  subj)));
  
  pow      = (abs(F*transpose(freq_shift.fourierspctrm)).^2)*P;
  pow      = log10(pow);
  s = std(pow, [], 2);
  u = mean(pow, 2);
  pow = (pow-repmat(u, [1 size(pow,2)]))./repmat(s, [1 size(pow,2)]);
  
  rho = corr(data_onset.trialinfo(:,7), pow', 'type', 'spearman'); %MvE
  rho=rho';
  
  filename = sprintf('/project/3011085.02/analysis/corr/sub-%03d/sub-%03d_corr_3Dgamma_rt.mat', subj, subj);
  save(filename, 'rho')
end
if dostat_pow_erf
  if ~exist('GA'); GA = input('send out single subject analyses (0), or continue to group analysis (1)?'); end
  if ~exist('whichFreq'); whichFreq = input('gamma (1), lowfreq (2)?'); end
  if ~GA
    sel = setdiff(1:33,10);
    for k = sel(:)'
      if whichFreq==1
        qsubfeval('erfosc_execute_pipeline','erfosc_script_jm',k,{'docorrpow_lcmv',1}, {'dofreq_short', 1},{'savefreq', 0},'memreq',8*1024^3,'timreq',59*60,'batchid',sprintf('subj%03d',k));
      elseif whichFreq==2
        qsubfeval('erfosc_execute_pipeline','erfosc_script_jm',k,{'docorrpow_lcmv_lowfreq',1}, {'dofreq_short_lowfreq', 1},{'savefreq', 0},'memreq',8*1024^3,'timreq',59*60,'batchid',sprintf('subj%03d',k));
      end
    end
  else
    load atlas_subparc374_8k.mat
    
    datadir = [project_dir 'analysis/'];
    erfosc_datainfo;
    k=1;
    for subj = allsubs
      if whichFreq==1
        filename = fullfile([datadir, 'corr/', sprintf('sub-%03d/sub-%03d_corrpowlcmv_gamma.mat', subj, subj)]);
      elseif whichFreq==2
        filename = fullfile([datadir, 'corr/', sprintf('sub-%03d/sub-%03d_corrpowlcmv_low.mat', subj, subj)]);
      end
      load(filename,'source');
      if k==1
        S=source;
        S.rho = source.rho;
        if whichFreq==1
          S_pupild.rho = source.partialrho(:,1);
          S_xy.rho = source.partialrho(:,2);
          S_eye.rho = source.partialrho(:,3);
        end
      else
        S.rho(:,k)=source.rho;
        if whichFreq==1
          S_pupild.rho(:,k) = source.partialrho(:,1);
          S_xy.rho(:,k) = source.partialrho(:,2);
          S_eye.rho(:,k) = source.partialrho(:,3);
        end
      end
      k=k+1;
    end
    S.rho = S.rho';
    exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
    S.rho(:, exclude_label) = nan; %MvE
    S.dimord = 'rpt_chan_freq';
    S.freq = 0;
    
    
    S2 = S;
    S2.rho(:) = 0;
    
    n = 32;
    
    cfgs = [];
    cfgs.method='montecarlo';
    cfgs.design=[ones(1,n) ones(1,n)*2;1:n 1:n];
    cfgs.statistic='ft_statfun_wilcoxon';
    cfgs.numrandomization=10000;
    cfgs.ivar=1;
    cfgs.uvar=2;
    cfgs.parameter='rho';
    cfgs.correctm='cluster';
    cfgs.clusterthreshold='nonparametric_individual';
    cfgs.connectivity = parcellation2connectivity_midline(atlas);
    cfgs.neighbours = cfgs.connectivity; %MvE
    cfgs.tail = 1;
    cfgs.clustertail = 1;
    cfgs.correcttail = 'prob';
    cfgs.clusteralpha = 0.05;
    cfgs.alpha = 0.05;
    
    stat=ft_freqstatistics(cfgs,S,S2);
    
    if whichFreq==1;
      S_pupild.rho = S_pupild.rho';
      S_xy.rho = S_xy.rho';
      S_eye.rho = s_eye.rho';
      S_pupild.rho(:, exclude_label) = nan;
      S_xy.rho(:, exclude_label) = nan;
      S_xy.rho(:, exclude_label) = nan;
      
      stat_pupild=ft_freqstatistics(cfgs,S_pupild,S2);
      stat_xy=ft_freqstatistics(cfgs,S_xy,S2);
      stat_eye=ft_freqstatistics(cfgs,S_eye,S2);
    end
    
    
    if whichFreq==1;
      filename = '/project/3011085.02/analysis/stat_peakpicking_gamma.mat';
      save(filename, 'stat', 'S', 'S_pupild', 'S_xy', 'stat_eye', 'stat_pupild', 'stat_xy');
    elseif whichFreq==2
      filename = '/project/3011085.02/analysis/stat_peakpicking_lowfreq.mat';
      save(filename, 'stat','S');
    end
  end
end
if dostat_erf_rt
  erfosc_datainfo;
  load atlas_subparc374_8k
  
  datadir = [project_dir 'analysis/'];
  k=1;
  
  d = dir('*corrpowlcmv_peakpicking_gamma.mat');
  
  
  for subj = allsubs
    filename = fullfile([datadir, 'corr/', sprintf('sub-%03d/sub-%03d_corrpowlcmv_gamma.mat', subj, subj)]);
    tmp = load(filename,'X', 'pow');
    X{k}=tmp.X;
    pow{k}=tmp.pow;
    k=k+1;
  end
  k=1;
  for subj=allsubs
    tmp = load(sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_rt.mat', subj,subj));
    rt{k} = tmp.rt;
    k=k+1;
  end
  
  for k=1:32
    rho{k} = corr(X{k}', rt{k}, 'type', 'spearman'); %MvE
  end
  
  exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
  
  source = [];
  source.brainordinate = atlas;
  source.label         = atlas.parcellationlabel;
  source.freq          = 0;
  source.dimord        = 'rpt_chan_freq';
  source.rho           = zeros(32,374);
  
  indx = 1:374;
  indx(exclude_label) = [];
  source.rho(:,indx) = cat(2,rho{:})';
  source.rho(:,exclude_label) = nan;
  
  
  ref=source;
  ref.rho(:) = 0;
  
  n = 32;
  
  cfgs = [];
  cfgs.method='montecarlo';
  cfgs.design=[ones(1,n) ones(1,n)*2;1:n 1:n];
  cfgs.statistic='ft_statfun_wilcoxon';
  cfgs.numrandomization=10000;
  cfgs.ivar=1;
  cfgs.uvar=2;
  cfgs.parameter='rho';
  cfgs.correctm='cluster';
  cfgs.clusterthreshold='nonparametric_individual';
  cfgs.connectivity = parcellation2connectivity_midline(atlas);
  cfgs.neighbours = cfgs.connectivity; %MvE
  cfgs.correcttail = 'prob';
  cfgs.tail = -1;
  cfgs.clustertail = -1;
  cfgs.alpha = 0.05;
  cfgs.clusteralpha = 0.01;
  
  stat=ft_freqstatistics(cfgs,source, ref);
  
  filename = '/project/3011085.02/analysis/stat_corr_peakpicking_rt.mat';
  save(filename, 'stat', 'source');
end