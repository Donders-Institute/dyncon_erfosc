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
    erf_osc_datainfo;
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
if ~exist('docorrelation', 'var'), docorrelation = false; end
if ~exist('docorrelation_lowfreq', 'var'), docorrelation_lowfreq = false; end
if ~exist('dolcmv_parc', 'var'),  dolcmv_parc = false; end
if ~exist('dolcmv_parc_msmall', 'var'),  dolcmv_parc_msmall = false; end
if ~exist('dolcmv_norm', 'var'),  dolcmv_norm = false; end
if ~exist('dolcmv_norm_msmall', 'var'),  dolcmv_norm_msmall = false; end
if ~exist('dosplitpow_lcmv', 'var'), dosplitpow_lcmv = false; end
if ~exist('docorrpow_lcmv', 'var'), docorrpow_lcmv = false; end
if ~exist('doPlanar', 'var'), doPlanar = false; end
if ~exist('doglm', 'var'), doglm = false; end
if ~exist('dosplitpow_source', 'var'), dosplitpow_source = false; end
if ~exist('doparcel_erf', 'var'), doparcel_erf = false; end
if ~exist('doresplocked', 'var'), doresplocked = false; end
if ~exist('docorrpow_lcmv_lowfreq', 'var'), docorrpow_lcmv_lowfreq = false; end
if ~exist('docorr_lcmv_eye', 'var'), docorr_lcmv_eye = false; end

if doparcel_erf, dolcmv_parc = true; end
if dodics_gamma, dofreq  = true; end
if dodics_lowfreq, dofreq  = true; end
if dofreq,       getdata = true; end
if docomp,       getdata = true; end
if dofreq_short, getdata = true; end
if dofreq_short_lowfreq, getdata = true; end
if dolcmv_parc,  getdata = true; end
if dolcmv_parc_msmall, getdata = true; end
if dolcmv_norm,  getdata = true; end
if dosplitpow_lcmv, getdata = true; end
if docorrpow_lcmv, getdata = true; end
if docorrpow_lcmv_lowfreq, getdata = true; end
if docorr_lcmv_eye, getdata = true; end
if doglm,     getdata = true; dolcmv_parc = true; end
if dosplitpow_source, getdata = true; end


% this chunk creates 2 data structures [-0.75 0.5]
if getdata
    if ~exist('undocomp','var')
        undocomp = false;
    end
    [p,f,e]       = fileparts(subject.dataset);
    basedir       = strrep(p, 'raw', 'processed');
    filename_data = fullfile(basedir, 'cleandata.mat');
    load(filename_data);
    
    if undocomp
        filename_comp = fullfile(basedir, 'icaComp.mat');
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

% this chunk does a dss decomposition on the onset-aligned data, and
% applies the topographies to the shift-aligned data
if docomp
    [comp_onset, comp_shift, comp_shift2] = erfosc_dss_onset(data_onset, data_shift);
    
    savedir = '/home/language/jansch/erfosc';
    save(fullfile(savedir, sprintf('sub-%03d_comp', subj)), 'comp_onset', 'comp_shift', 'comp_shift2');
end

% this chunk does spectral decomposition
if dofreq_short
    peakpicking;%MvE
    latency = peaks(subj,1).*[1 1] - 0.02 - [0.20 1./600];%MvE
    
    foi     = subject.gammapeak(end).*[1 1];
    smo     = max(10, diff(subject.gammaband(end,:))./2);
    [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    if ~exist('savefreq', 'var')
        savefreq = false;
    end
    if savefreq
        savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';%MvE
        save(fullfile(savedir, sprintf('sub-%03d_freqshort_mve', subj)), 'freq_shift', 'P', 'latency');
    end
end

% this chunk does spectral decomposition
if dofreq_short_lowfreq
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow_low.mat',subj));
    peakpicking;
    latency = peaks(subj,1).*[1 1] - 0.02 - [0.4 1./600];%MvE
    foi = [peakfreq peakfreq];
    smo     = 2.5;
    [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    if ~exist('savefreq', 'var')
        savefreq = false;
    end
    if savefreq
        savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
        save(fullfile(savedir, sprintf('sub-%03d_freqshort_low', subj)), 'freq_shift', 'P', 'latency');
    end
end

% this chunk does spectral decomposition
if dofreq
    
    if ~exist('latency', 'var')
        latency = [-inf 0-1./data_onset.fsample];
    end
    if dodics_lowfreq
        if ~exist('peakfreq', 'var')
            load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow_low.mat',subj));
        end
        foi = [peakfreq peakfreq]; smo = 2;
        latency = [-0.75+1./data_onset.fsample 0-1./data_onset.fsample];
        [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    else
        [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject);
    end
    
    if ~exist('savefreq', 'var')
        savefreq = false;
    end
    if savefreq
        %         savedir = '/home/language/jansch/erfosc';
        savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
        save(fullfile(savedir, sprintf('sub-%03d_freq', subj)), 'freq_shift', 'P', 'latency');
    end
    
end

% this chunk does source reconstruction
if dodics_gamma
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    [source_onset, source_shift, Tval, F] = erfosc_dics_gamma(freq_onset, freq_shift, headmodel, sourcemodel);
    source_onset = rmfield(source_onset, 'cfg');
    source_shift = rmfield(source_shift, 'cfg');
    
    savedir = '/home/language/jansch/erfosc';
    save(fullfile(savedir, sprintf('sub-%03d_source', subj)), 'source_onset', 'source_shift', 'Tval', 'F');
end

if dodics_lowfreq
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    [source_onset, source_shift, Tval, F] = erfosc_dics_alpha(freq_onset, freq_shift, headmodel, sourcemodel);
    source_onset = rmfield(source_onset, 'cfg');
    source_shift = rmfield(source_shift, 'cfg');
    
    savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
    save(fullfile(savedir, sprintf('sub-%03d_source_low', subj)), 'source_shift', 'Tval', 'F');
end

if dolcmv_parc
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    load('atlas_subparc374_8k.mat');
    if doresplocked
        [source_parc] = erfosc_lcmv_parc(data_resp, headmodel, sourcemodel, atlas, doresplocked);
    else
        [source_parc] = erfosc_lcmv_parc(data_shift, headmodel, sourcemodel, atlas, doresplocked);
        savedir = '/home/language/jansch/erfosc';
        save(fullfile(savedir, sprintf('sub-%03d_lcmv', subj)), 'source_parc');
    end
end

if dolcmv_norm
    datadir  = '/home/language/jansch/erfosc';
    filename = fullfile(datadir, sprintf('sub-%03d_lcmv', subj));
    load(filename);
    
    tmpcfg = [];
    tmpcfg.latency = [-0.2 -1./600];
    data_shift = ft_selectdata(tmpcfg, data_shift);
    
    tmpcfg = [];
    tmpcfg.covariance = 'yes';
    tmp = ft_timelockanalysis(tmpcfg, data_shift);
    
    noise = zeros(numel(source_parc.label),1);
    for k = 1:numel(source_parc.label)
        tmpF = source_parc.F{k}(1,:);
        tmpC = sqrt(tmpF*tmp.cov*tmpF');
        noise(k) = tmpC;
    end
    save(filename, 'noise', '-append');
    
end

if dolcmv_parc_msmall
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    load('atlas_MSMAll_8k_subparc.mat');
    [source_parc] = erfosc_lcmv_parc(data_shift, headmodel, sourcemodel, atlas);
    
    savedir = '/home/language/jansch/erfosc';
    save(fullfile(savedir, sprintf('sub-%03d_lcmv_msmall', subj)), 'source_parc');
end

if dolcmv_norm_msmall
    datadir  = '/home/language/jansch/erfosc';
    filename = fullfile(datadir, sprintf('sub-%03d_lcmv_msmall', subj));
    load(filename);
    
    tmpcfg = [];
    tmpcfg.latency = [-0.2 -1./600];
    data_shift = ft_selectdata(tmpcfg, data_shift);
    
    tmpcfg = [];
    tmpcfg.covariance = 'yes';
    tmp = ft_timelockanalysis(tmpcfg, data_shift);
    
    noise = zeros(numel(source_parc.label),1);
    for k = 1:numel(source_parc.label)
        tmpF = source_parc.F{k}(1,:);
        tmpC = sqrt(tmpF*tmp.cov*tmpF');
        noise(k) = tmpC;
    end
    save(filename, 'noise', '-append');
    
end

% this chunk extracts single trial power and amplitude of
if docorrelation
    %     erfosc_comppeaks;
    erfosc_lcmvpeaks
    datadir1 = '/home/language/jansch/erfosc';
    datadir2 = '/project/3011085.02/scripts/erfosc/analysis_JM_data/'
    load(fullfile(datadir1, sprintf('sub-%03d_source',    subj)));
    %     load(fullfile(datadir, sprintf('sub-%03d_comp',      subj)));
    load(fullfile(datadir2, sprintf('sub-%03d_splitpow', subj)), 'datapst')
    load(fullfile(datadir1, sprintf('sub-%03d_lcmv', subj)));
    load(fullfile(datadir2, sprintf('sub-%03d_freqshort', subj)));
    
    [val1 idx1] = sort(datapst.trialinfo(:,1));
    datapst.trial = datapst.trial(idx1);
    datapst.trialinfo = datapst.trialinfo(idx1,:);
    
    [m, idx] = max(Tval);
    pow      = (abs(F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    %pow      = log10(pow);
    
    %     erf      = cellrowselect(comp_shift.trial, comp_id(subj));
    %     erf      = cat(1,erf{:});
    %     erf      = polarity(subj).*ft_preproc_baselinecorrect(erf, 391, 451);
    erf = cat(1, datapst.trial{:});
    
    %     ix1 = nearest(comp_shift.time{1}, peaks(subj,1));
    %     ix2 = nearest(comp_shift.time{1}, peaks(subj,2));
    ix1 = nearest(datapst.time{1}, peaks(subj,1));
    ix2 = nearest(datapst.time{1}, peaks(subj,2));
    signpeak = sign(mean(mean(erf(:,ix1:ix2),2)));
    erf = signpeak.*erf;
    amp = mean(erf(:,ix1:ix2), 2);
    
    [rho, pval] = corr(log10(pow(:)), amp(:), 'type', 'spearman');
    
    save(fullfile(datadir2, sprintf('sub-%03d_corr', subj)), 'amp', 'pow', 'rho', 'pval', 'erf');
end

if docorrelation_lowfreq
    erfosc_lcmvpeaks;
    datadir1 = '/home/language/jansch/erfosc';
    datadir2 = '/project/3011085.02/scripts/erfosc/analysis_JM_data/'
    load(fullfile(datadir1, sprintf('sub-%03d_source_alpha',    subj)));
    load(fullfile(datadir2, sprintf('sub-%03d_splitpow', subj)), 'datapst')
    load(fullfile(datadir1, sprintf('sub-%03d_lcmv', subj)));
    load(fullfile(datadir1, sprintf('sub-%03d_freqshort_alpha', subj)));
    
    [val1 idx1] = sort(datapst.trialinfo(:,1));
    datapst.trial = datapst.trial(idx1);
    datapst.trialinfo = datapst.trialinfo(idx1,:);
    
    [m, idx] = min(Tval);
    pow      = (abs(F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;

    erf = cat(1, datapst.trial{:});
    
    ix1 = nearest(datapst.time{1}, peaks(subj,1));
    ix2 = nearest(datapst.time{1}, peaks(subj,2));
    signpeak = sign(mean(mean(erf(:,ix1:ix2),2)));
    erf = signpeak.*erf;
    

    amp = mean(erf(:,ix1:ix2), 2);
    
    [rho, pval] = corr(log10(pow(:)), amp(:), 'type', 'spearman');
    
    save(fullfile(datadir2, sprintf('sub-%03d_corr_alpha', subj)), 'amp', 'pow', 'rho', 'pval', 'erf');
end

if dosplitpow_lcmv
    erfosc_lcmvpeaks;
    
    datadir = '/home/language/jansch/erfosc';
    load(fullfile(datadir, sprintf('sub-%03d_lcmv',    subj)));
    source_parc.avg = diag(1./noise)*source_parc.avg;
    
    ix1 = nearest(source_parc.time, peaks(subj,1));
    ix2 = nearest(source_parc.time, peaks(subj,2));
    
    tmpcfg = [];
    tmpcfg.latency = [-0.1 0.5-1./600];
    datapst = ft_selectdata(tmpcfg, data_shift);
    tmpcfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02; %0.01;
    datapre = ft_selectdata(tmpcfg, data_shift);
    %clear data_shift;
    
    source_parc.avg  = ft_preproc_baselinecorrect(source_parc.avg, 1, 60);
    [maxval, maxidx] = max(abs(mean(source_parc.avg(:,ix1:ix2),2)));
    signpeak         = sign(mean(source_parc.avg(maxidx,ix1:ix2),2));
    fprintf('parcel with max amplitude = %s\n', source_parc.label{maxidx});
    
    plot_alternative_parcel(source_parc.label{maxidx});
    datapst.trial = source_parc.F{maxidx}(1,:)*datapst.trial;
    datapst.label = source_parc.label(maxidx);
    
    tmpcfg = [];
    tmpcfg.demean = 'yes';
    tmpcfg.baselinewindow = [-inf 0];
    datapst = ft_preprocessing(tmpcfg, datapst);
    
    tmp        = cat(1, datapst.trial{:}).*signpeak; % meak peak up-going
    tmp        = ft_preproc_baselinecorrect(tmp, 1, 60);
    [srt, idx] = sort(mean(tmp(:,ix1:ix2),2));
    datapst.trial = datapst.trial(idx);
    datapst.trialinfo = datapst.trialinfo(idx,:);
    datapre.trial = datapre.trial(idx);
    datapre.trialinfo = datapre.trialinfo(idx,:);
    data_shift.trial = data_shift.trial(idx);
    data_shift.trialinfo = data_shift.trialinfo(idx,:);
    
    
    % powerspectrum
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.tapsmofrq = 10;
    cfg.foilim = [0 120];
    cfg.trials = 1:100;
    cfg.pad = 1;
    
    if doPlanar
        cfgb                 = [];
        cfgb.method          = 'template';
        cfgb.template        = 'CTF275_neighb.mat';
        cfgb.neighbours      = ft_prepare_neighbours(cfgb, datapre);
        cfgb.method          = 'sincos';
        datapre_planar         = ft_megplanar(cfgb, datapre);
        cfgb = [];
        cfgb.method = 'sum';
        f1 = ft_combineplanar(cfgb, ft_freqanalysis(cfg, datapre_planar));
        cfg.trials = numel(datapre.trial)-100 + (1:100);
        f2 = ft_combineplanar(cfgb, ft_freqanalysis(cfg, datapre_planar));
    else
        f1 = ft_freqanalysis(cfg, datapre);
        cfg.trials = numel(datapre.trial)-100 + (1:100);
        f2 = ft_freqanalysis(cfg, datapre);
    end
    
    % TFR
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.output = 'pow';
    cfg.pad    = 2;
    cfg.foi    = [2:2:120];
    cfg.tapsmofrq = ones(1,numel(cfg.foi)).*8;
    cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.25;
    cfg.toi    = (-150:6:90)./300;
    cfg.trials = 1:100;
    
    if doPlanar
        cfgb                 = [];
        cfgb.method          = 'template';
        cfgb.template        = 'CTF275_neighb.mat';
        cfgb.neighbours      = ft_prepare_neighbours(cfgb, data_shift);
        cfgb.method          = 'sincos';
        data_shift_planar         = ft_megplanar(cfgb, data_shift);
        cfgb = [];
        cfgb.method = 'sum';
        tfr1 = ft_combineplanar(cfgb, ft_freqanalysis(cfg, data_shift_planar));
        cfg.trials = numel(data_shift.trial)-100 + (1:100);
        tfr2 = ft_combineplanar(cfgb, ft_freqanalysis(cfg, data_shift_planar));
    else
        tfr1 = ft_freqanalysis(cfg, data_shift);
        cfg.trials = numel(data_shift.trial)-100 + (1:100);
        tfr2 = ft_freqanalysis(cfg, data_shift);
    end
    
    datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
    if doPlanar
        save(fullfile(datadir, sprintf('sub-%03d_splitpow_plcmb', subj)), 'f1', 'f2', 'datapst','idx');
        save(fullfile(datadir, sprintf('sub-%03d_splitpow_tfr_plcmb', subj)), 'tfr1', 'tfr2');
    else
        save(fullfile(datadir, sprintf('sub-%03d_splitpow', subj)), 'f1', 'f2', 'datapst');
        save(fullfile(datadir, sprintf('sub-%03d_splitpow_tfr', subj)), 'tfr1', 'tfr2');
    end
end
if doglm
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaPow');
    
    if doresplocked
        parceldata_shift = data_resp;
    else
        parceldata_shift = data_shift;
    end
    
    parceldata_shift.label = source_parc.label
    for k=1:length(source_parc.F);
        spatfilter_parcel(k,:) = source_parc.F{k}(1,:);
    end
    parceldata_shift.trial = spatfilter_parcel*parceldata_shift.trial;
    
    cfg=[];
    cfg.lpfilter = 'yes';
    cfg.lpfilttype = 'firws';
    cfg.lpfreq = 30;
    cfg.lpfiltdir = 'onepass-reverse-zerophase';
    cfg.preproc.demean = 'yes';
    if ~doresplocked
        cfg.preproc.baselinewindow = [-0.1 0];
    end
    parceldata_shift = ft_preprocessing(cfg, parceldata_shift);
    
    if doresplocked
        % add nans for short trials
        tend = zeros(numel(parceldata_shift.time),1);
        for k=1:numel(parceldata_shift.time);
            tend(k,1) = parceldata_shift.time{k}(end); % save the last time point for each trial
            parceldata_shift.time{k} = [parceldata_shift.time{k}, parceldata_shift.time{k}(end)+1/parceldata_shift.fsample:1/parceldata_shift.fsample:0.4];
            parceldata_shift.trial{k} = [parceldata_shift.trial{k}, nan(numel(parceldata_shift.label), numel(parceldata_shift.time{k})-size(parceldata_shift.trial{k},2))];
        end
        cfg=[];
        cfg.latency = [-0.5 0.4];
        erfdata = ft_selectdata(cfg, parceldata_shift);
        tlck=rmfield(erfdata,{'trial', 'time', 'trialinfo'});
        tlck.avg = nanmean(cat(3,erfdata.trial{:}),3);
        tlck.dimord = 'chan_time';
        tlck.time = erfdata.time{1};
    else
        cfg=[];
        cfg.latency = [-0.1 0.5];
        erfdata = ft_selectdata(cfg, parceldata_shift);
        
        cfg=[];
        cfg.latency = [-0.1 0];
        tmp = ft_selectdata(cfg, erfdata);
        tmp.trial = cat(2, tmp.trial{:});
        baseline_std = std(tmp.trial, [], 2);
        cfg=[];
        cfg.vartrllength = 2;
        tlck = ft_timelockanalysis(cfg, erfdata);
    end
    
    
    
    erfdata.trial = cat(3,erfdata.trial{:});
    erfdata.trial = permute(erfdata.trial, [3,1,2]);
    erfdata.time = erfdata.time{1};
    
    design = [gammaPow;((data_shift.trialinfo(:,5)-data_shift.trialinfo(:,4))/1200)'; ones(size(gammaPow))];
    cfg=[];
    cfg.glm.statistic = 'beta';
    cfg.glm.standardise = false;
    
    for k=1:length(erfdata.label)
        dat = [squeeze(erfdata.trial(:,k,:))]';
        if doresplocked
            dat = (dat - repmat(nanmean(dat,2),[1 length(erfdata.trialinfo)]));
            dat(isnan(dat))=0;
        else
            dat = (dat - repmat(mean(dat,2),[1 length(erfdata.trialinfo)]));
        end
        tmp = statfun_glm(cfg, dat, design);
        betas_tmp(k,:) = tmp.stat(:,1);
    end
    
    betas        = rmfield(erfdata,{'trial', 'cfg'});
    betas.avg    = betas_tmp;
    betas.time   = erfdata.time;
    betas.dimord = 'chan_time';
    
    savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
    if doresplocked
        save(fullfile([savedir, sprintf('sub-%03d_glm_parcelresp', subj)]), 'betas', 'tlck','tend');
    else
        save(fullfile([savedir, sprintf('sub-%03d_parcel_blstd', subj)]), 'baseline_std');
        save(fullfile([savedir, sprintf('sub-%03d_glm_parcel', subj)]), 'betas', 'tlck');
    end
end
if dosplitpow_source
    datadir1 = '/home/language/jansch/erfosc';
    datadir2 = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
    load(fullfile(datadir2, sprintf('sub-%03d_splitpow_plcmb', subj)), 'datapst'); % take the trial indexes of sorted trials
    
    load(fullfile(datadir1, 'atlas_subparc374_8k.mat'))
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    foi = subject.gammaband;
    [source_parc] = erfosc_lcmv_parc_gamma(data_shift, headmodel, sourcemodel, atlas, foi);
    
    parceldata_shift = data_shift;
    
    parceldata_shift.label = source_parc.label
    for k=1:length(source_parc.F);
        spatfilter_parcel(k,:) = source_parc.F{k}(1,:);
    end
    parceldata_shift.trial = spatfilter_parcel*parceldata_shift.trial;
    
    % sort based on erf amplitude (as in datapst)
    erftrlinfo = datapst.trialinfo(:,1);
    
    parceltrlinfo=parceldata_shift.trialinfo(:,1);
    [~, order] = ismember(erftrlinfo, parceltrlinfo);
    parceldata_shift.trialinfo = parceldata_shift.trialinfo(order, :); % sort parceldata on erf amp (low-high)
    parceldata_shift.trial = parceldata_shift.trial(order);
    
    % TFR
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.output = 'pow';
    cfg.pad    = 2;
    cfg.foi    = 20:2:100
    cfg.tapsmofrq = ones(1,numel(cfg.foi)).*8;
    cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.25;
    cfg.toi    = (-300:6:150)./300;
    cfg.trials = 1:100;
    tfr_high_1 = ft_freqanalysis(cfg, parceldata_shift); % low ERF amp
    cfg.trials = numel(data_shift.trial)-100 + (1:100);
    tfr_high_2 = ft_freqanalysis(cfg, parceldata_shift); % high ERF amp
    
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.output = 'pow';
    cfg.pad    = 2;
    cfg.foi    = 2:2:30
    cfg.taper  = 'hanning';
    cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.25;
    cfg.toi    = (-300:6:150)./300;
    cfg.trials = 1:100;
    tfr_low_1 = ft_freqanalysis(cfg, parceldata_shift); % low ERF amp
    cfg.trials = numel(data_shift.trial)-100 + (1:100);
    tfr_low_2 = ft_freqanalysis(cfg, parceldata_shift); % high ERF amp
    
    savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
    save(fullfile(savedir, sprintf('sub-%03d_splitpow_source', subj)), 'tfr_low_1','tfr_low_2','tfr_high_1', 'tfr_high_2')
end

if doparcel_erf
    parceldata_shift = data_shift;
    
    parceldata_shift.label = source_parc.label;
    for k=1:length(source_parc.F);
        spatfilter_parcel(k,:) = source_parc.F{k}(1,:);
    end
    parceldata_shift.trial = spatfilter_parcel*parceldata_shift.trial;
    
    [~, idx] = sort(parceldata_shift.trialinfo(:,6)-parceldata_shift.trialinfo(:,5), 'ascend');
    parceldata_shift.trial = parceldata_shift.trial(idx);
    parceldata_shift.trialinfo = parceldata_shift.trialinfo(idx);
    
    savedir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
    save(fullfile(savedir, sprintf('sub-%03d_erfparc', subj)), 'parceldata_shift');
    
end


if docorrpow_lcmv
    peakpicking;
    
    datadir = '/home/language/jansch/erfosc';
    load(fullfile(datadir, sprintf('sub-%03d_lcmv',    subj)));
    source_parc.avg = diag(1./noise)*source_parc.avg;
    
    ix1 = nearest(source_parc.time, peaks(subj,1));
    ix2 = nearest(source_parc.time, peaks(subj,2));
    
    tmpcfg = [];
    tmpcfg.latency = [-0.1 0.5-1./600];
    datapst = ft_selectdata(tmpcfg, data_shift);
    tmpcfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02; 
    datapre = ft_selectdata(tmpcfg, data_shift);
    
    source_parc.avg  = ft_preproc_baselinecorrect(source_parc.avg, 1, 60);
    [maxval, maxidx] = max(abs(mean(source_parc.avg(:,ix1:ix2),2)));
    signpeak         = sign(mean(source_parc.avg(maxidx,ix1:ix2),2));
    fprintf('parcel with max amplitude = %s\n', source_parc.label{maxidx});
    
    for k = 1:numel(source_parc.label)
        F(k,:) = source_parc.F{k}(1,:);
    end
    datapst.trial = F*datapst.trial;
    datapst.label = source_parc.label;
    
    tmpcfg = [];
    tmpcfg.demean = 'yes';
    tmpcfg.baselinewindow = [-inf 0];
    tmpcfg.lpfilter = 'yes';
    tmpcfg.lpfreq = 30;
    tmpcfg.lpfilttype = 'firws';
    tmpcfg.lpfiltdir = 'onepass-reverse-zerophase'
    datapst = ft_preprocessing(tmpcfg, datapst);
    
    tmpcfg = [];
    tmpcfg.latency = peaks(subj,:);%JM
    tmpcfg.avgovertime = 'yes';
    datapeak = ft_selectdata(tmpcfg,datapst);
    
    X = cat(2,datapeak.trial{:});
    Xpow = abs(mean(X,2));
    signswap = diag(sign(mean(X,2)));
    X = signswap*X; % let the amplitude be on average positive
    X = standardise(X,2);
    
    tmpcfg = [];
    tlckpst = ft_timelockanalysis(tmpcfg, datapst);
    tlckpst.avg = signswap*tlckpst.avg;
    
    load(fullfile(datadir, sprintf('sub-%03d_source',    subj)));
%     load(fullfile(datadir, sprintf('sub-%03d_freqshort',    subj)));
    [m, idx] = max(Tval);
    pow      = (abs(F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    pow      = standardise(log10(pow(:)));
    
    rho = corr(X', pow, 'type', 'spearman'); %MvE

    tmp = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/eyedata.mat', subj));
    % get gamma-ERF correlation, accounting for pupil diameter, without confounding eye position.
    cfg=[];
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
    eye = ft_timelockanalysis(cfg,tmp.data_shift);
    cfg=[];
    cfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    pupild = erfosc_regress_eye(ft_selectdata_new(cfg, eye), {'UADC007'}, {'visAngleX', 'visAngleY'});
    cfg=[];
    cfg.avgovertime = 'yes';
    pupild = ft_selectdata(cfg, pupild);
    pupild = standardise(pupild.trial);
    partialrho1 = partialcorr(X', pow, pupild, 'type', 'spearman'); 
    
    % get correlation gamma-ERF, accounting for eye position, without confound pupil diameter
    idx = match_str(eye.label, {'visAngleX', 'visAngleY'});
    eye.trial(:,end+1,:) = (eye.trial(:,idx(1),:).^2 + eye.trial(:,idx(2),:).^2).^0.5;
    eye.label{end+1} = 'distance';
    cfg=[];
    cfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    eyepos = erfosc_regress_eye(ft_selectdata_new(cfg, eye), {'distance'}, {'UADC007'});
    cfg=[];
    cfg.avgovertime = 'yes';
    distance = ft_selectdata(cfg, eyepos);
    distance = standardise(distance.trial);    
    partialrho2 = partialcorr(X', pow, distance, 'type', 'spearman');
    
    % get correlation gamma-ERF, accounting for eye position and pupil
    % diameter (both not confounded by the other)
    partialrho3 = partialcorr(X', pow, [distance, pupild], 'type', 'spearman'); 
        
    load atlas_subparc374_8k
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
    source.rho(indx)    = rho;
    source.partialrho(indx,1) = partialrho1;
    source.partialrho(indx,2) = partialrho2;
    source.partialrho(indx,3) = partialrho3;
    source.pow(indx)    = Xpow(:);
    
    datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
    save(fullfile(datadir, sprintf('sub-%03d_corrpowlcmv_gamma', subj)), 'source', 'pow', 'X', 'tlckpst', 'pupild', 'distance');
end
if docorrpow_lcmv_lowfreq
    peakpicking;
    
    datadir = '/home/language/jansch/erfosc';
    load(fullfile(datadir, sprintf('sub-%03d_lcmv',    subj)));
    source_parc.avg = diag(1./noise)*source_parc.avg;
    
    ix1 = nearest(source_parc.time, peaks(subj,1));
    ix2 = nearest(source_parc.time, peaks(subj,2));
    
    tmpcfg = [];
    tmpcfg.latency = [-0.1 0.5-1./600];
    datapst = ft_selectdata(tmpcfg, data_shift);
    tmpcfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    
    source_parc.avg  = ft_preproc_baselinecorrect(source_parc.avg, 1, 60);
    [maxval, maxidx] = max(abs(mean(source_parc.avg(:,ix1:ix2),2)));
    signpeak         = sign(mean(source_parc.avg(maxidx,ix1:ix2),2));
    fprintf('parcel with max amplitude = %s\n', source_parc.label{maxidx});
    
    for k = 1:numel(source_parc.label)
        F(k,:) = source_parc.F{k}(1,:);
    end
    datapst.trial = F*datapst.trial;
    datapst.label = source_parc.label;
    
    tmpcfg = [];
    tmpcfg.demean = 'yes';
    tmpcfg.baselinewindow = [-inf 0];
    tmpcfg.lpfilter = 'yes';
    tmpcfg.lpfreq = 30;
    tmpcfg.lpfilttype = 'firws';
    tmpcfg.lpfiltdir = 'onepass-reverse-zerophase';
    datapst = ft_preprocessing(tmpcfg, datapst);
    
    tmpcfg = [];
    tmpcfg.latency = peaks(subj,:);
    tmpcfg.avgovertime = 'yes';
    datapeak = ft_selectdata(tmpcfg,datapst);
    
    X = cat(2,datapeak.trial{:});
    Xpow = abs(mean(X,2));
    signswap = diag(sign(mean(X,2)));
    X = signswap*X; % let the amplitude be on average positive
    X = standardise(X,2);
    
    tmpcfg = [];
    tlckpst = ft_timelockanalysis(tmpcfg, datapst);
    tlckpst.avg = signswap*tlckpst.avg;
    
    datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
    load(fullfile(datadir, sprintf('sub-%03d_source_low',    subj)));
    load(fullfile(datadir, sprintf('sub-%03d_freqshort_low', subj)));
    [m, idx] = min(Tval);
    pow      = (abs(F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    pow      = standardise(log10(pow(:)));
    
    rho = corr(X', pow, 'type', 'spearman'); %MvE
    load atlas_subparc374_8k
    exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
    
    source = [];
    source.brainordinate = atlas;
    source.label         = atlas.parcellationlabel;
    source.rho           = zeros(374,1);
    source.pow           = zeros(374,1);
    source.dimord        = 'chan';
    
    indx = 1:374;
    indx(exclude_label) = [];
    source.rho(indx)    = rho;
    source.pow(indx)    = Xpow(:);
    
    
    datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
    save(fullfile(datadir, sprintf('sub-%03d_corrpowlcmv_low', subj)), 'source', 'pow', 'X', 'tlckpst');
end
if docorr_lcmv_eye
    peakpicking;
    
    datadir = '/home/language/jansch/erfosc';
    load(fullfile(datadir, sprintf('sub-%03d_lcmv',    subj)));
    source_parc.avg = diag(1./noise)*source_parc.avg;
    
    ix1 = nearest(source_parc.time, peaks(subj,1));
    ix2 = nearest(source_parc.time, peaks(subj,2));
    
    tmpcfg = [];
    tmpcfg.latency = [-0.1 0.5-1./600];
    datapst = ft_selectdata(tmpcfg, data_shift);
    tmpcfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    
    source_parc.avg  = ft_preproc_baselinecorrect(source_parc.avg, 1, 60);
    [maxval, maxidx] = max(abs(mean(source_parc.avg(:,ix1:ix2),2)));
    signpeak         = sign(mean(source_parc.avg(maxidx,ix1:ix2),2));
    fprintf('parcel with max amplitude = %s\n', source_parc.label{maxidx});
    
    for k = 1:numel(source_parc.label)
        F(k,:) = source_parc.F{k}(1,:);
    end
    datapst.trial = F*datapst.trial;
    datapst.label = source_parc.label;
    
    tmpcfg = [];
    tmpcfg.demean = 'yes';
    tmpcfg.baselinewindow = [-inf 0];
    tmpcfg.lpfilter = 'yes';
    tmpcfg.lpfreq = 30;
    tmpcfg.lpfilttype = 'firws';
    tmpcfg.lpfiltdir = 'onepass-reverse-zerophase';
    datapst = ft_preprocessing(tmpcfg, datapst);
    
    tmpcfg = [];
    tmpcfg.latency = peaks(subj,:);
    tmpcfg.avgovertime = 'yes';
    datapeak = ft_selectdata(tmpcfg,datapst);
    
    X = cat(2,datapeak.trial{:});
    Xpow = abs(mean(X,2));
    signswap = diag(sign(mean(X,2)));
    X = signswap*X; % let the amplitude be on average positive
    X = standardise(X,2);
    
    tmpcfg = [];
    tlckpst = ft_timelockanalysis(tmpcfg, datapst);
    tlckpst.avg = signswap*tlckpst.avg;
    
    eye = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/eyedata.mat', subj), 'data_shift');
    eye = eye.data_shift;
    cfg=[];
    cfg.latency = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    eye = erfosc_regress_eye(ft_selectdata(cfg, eye), {'UADC007'}, {'visAngleX', 'visAngleY'});
    cfg=[];
    cfg.avgovertime = 'yes';
    eye = ft_selectdata(cfg, eye);
    pupild = eye.trial; 
    
    rho = corr(X', pupild, 'type', 'spearman');
    load atlas_subparc374_8k
    exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'}); %MvE
    
    source = [];
    source.brainordinate = atlas;
    source.label         = atlas.parcellationlabel;
    source.rho           = zeros(374,1);
    source.diameter      = zeros(374,1);
    source.dimord        = 'chan';
    
    indx = 1:374;
    indx(exclude_label) = [];
    source.rho(indx,1)  = rho;

    
    datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
    save(fullfile(datadir, sprintf('sub-%03d_corrlcmv_eye', subj)), 'source', 'rho','pupild');
end