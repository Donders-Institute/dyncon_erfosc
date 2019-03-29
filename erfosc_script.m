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

% this chunk does spectral decomposition
if dofreq_short || dofreq_short_lowfreq
    peakpicking;%MvE
    if dofreq_short
        latency = peaks(subj,1).*[1 1] - 0.02 - [0.20 1./600];%MvE
        foi     = subject.gammapeak(end).*[1 1];
        smo     = max(10, diff(subject.gammaband(end,:))./2);
    elseif dofreq_short_lowfreq
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
        if dofreq_short_lowfreq
            filename = [filename, '_low'];
        end
        save(filename, 'freq_shift', 'P', 'latency');
    end
end


% this chunk does spectral decomposition
if dofreq
    if ~exist('latency', 'var')
        latency = [-inf 0-1./data_onset.fsample];
    end
    if dodics_lowfreq
        if ~exist('peakfreq', 'var')
            load([project_dir, sprintf('analysis/freq/sub-%03d/sub-%03d_pow_low.mat',subj,subj)]);
        end
        foi = [peakfreq peakfreq]; smo = 2;
        latency = [-0.75+1./data_onset.fsample 0-1./data_onset.fsample];
        [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo);
    else
        [freq_onset, freq_shift, P, pow_onset, pow_shift] = erfosc_freq(data_onset, data_shift, latency, subject);
    end
    
    if ~exist('savefreq', 'var')
        savefreq = false;
    end
    if savefreq
        savedir = [project_dir 'analysis/freq'];
        save(fullfile(savedir, sprintf('sub-%03d/sub-%03d_freq', subj,subj)), 'freq_shift', 'P', 'latency');
    end
end

% this chunk does source reconstruction
if dodics_gamma || dodics_lowfreq
    load(fullfile(subject.mridir,'preproc','headmodel.mat'));
    cfg=[];
    cfg.comment = 'load subject`s headmodel';
    headmodel = ft_annotate(cfg, headmodel);
    load(fullfile(subject.mridir,'preproc','sourcemodel2d.mat'));
    cfg.comment = 'load subject`s 2D headmodel';
    sourcemodel = ft_annotate(cfg, sourcemodel);
    [source_onset, source_shift, Tval, F] = erfosc_dics(freq_onset, freq_shift, headmodel, sourcemodel);
    source_onset = rmfield(source_onset, 'cfg');
    source_shift = rmfield(source_shift, 'cfg');
    
    savedir = [project_dir 'analysis/source/'];
    if dosave
        filename = fullfile(savedir, sprintf('sub-%03d/sub-%03d_source', subj,subj));
        if dodics_lowfreq
            filename = [filename, '_low'];
        end
        save(filename, 'source_onset', 'source_shift', 'Tval', 'F');
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
%     datadir  = [project_dir 'analysis/source/'];
%     filename = fullfile(datadir, sprintf('sub-%03d/sub-%03d_lcmv', subj,subj));
%     load(filename);
    
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
    
    datadir = [project_dir 'analysis/'];
    load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_lcmv',    subj, subj)));
    cfg=[];
    cfg.matrix = repmat(1./source_parc.noise, [1 size(source_parc.avg,2)]);
    cfg.parameter = 'avg';
    cfg.operation = 'multiply';
    cfg.comment = 'devide the average by the noise';
    source_parc = ft_math(cfg, source_parc);
    
    cfg=[];
    cfg.baseline = [-0.1 -1./data_shift.fsample];
    source_parc = ft_timelockbaseline(cfg, source_parc);
    
    cfg = [];
    cfg.latency = [-0.1 0.5-1./data_shift.fsample];
    datapst = ft_selectdata(cfg, data_shift);
    
    datapst.trial = source_parc.F*datapst.trial;
    datapst.label = source_parc.label;
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

    cfg=[];
    cfg.operation = 'multiply';
    cfg.matrix = 
    
    X = datapeak.trial;
    Xpow = abs(mean(X,1));
    signswap = diag(sign(mean(X,1)));
    X = X*signswap; % let the amplitude be on average positive
    datapeak.amp = standardise(X,1);
    
    cfg=[];
    cfg.comment = 'Create .amp by aligning the parcels such that the trial-average of .trial is positive. Standardise over trials.';
    datapeak = ft_annotate(cfg, datapeak);
    
    cfg = [];
    datapst = ft_timelockanalysis(cfg, datapst);
    cfg = [];
    cfg.matrix = repmat(diag(signswap), [1 size(datapst.avg,2)]);
    cfg.operation = 'multiply';
    cfg.comment = sprintf('align the parcels such that the trial-average of the [%d %d]s window is positive.',peaks(subj,:));
    cfg.parameter = 'avg';
    datapst = ft_math(cfg, datapst);
%FIXME continue here.    
    load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_source',   subj, subj)));
    if ~exist('freq_shift')
        load(fullfile(datadir, 'freq/', sprintf('sub-%03d/sub-%03d_freqshort',  subj,  subj)));
    end
    [m, idx] = max(Tval);
    pow      = (abs(F(idx,:)*transpose(freq_shift.fourierspctrm)).^2)*P;
    pow      = standardise(log10(pow(:)));
    
    rho = corr(datapeak.amp, pow, 'type', 'spearman'); %MvE
    
    
    tmp = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_eyedata.mat', subj, subj));
    % get gamma-ERF correlation, accounting for pupil diameter, without confounding eye position.
    cfg=[];
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
    eye = ft_timelockanalysis(cfg,tmp.data_shift);
    cfg=[];
    cfg.toilim = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    pupild = erfosc_regress_eye(ft_redefinetrial(cfg, eye), {'UADC007'}, {'visAngleX', 'visAngleY'});
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
    cfg.toilim = [-0.2 -1./600] + peaks(subj,1) - 0.02;
    eyepos = erfosc_regress_eye(ft_redefinetrial(cfg, eye), {'distance'}, {'UADC007'});
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
    
    datadir = [project_dir 'analysis/corr/'];
    save(fullfile(datadir, sprintf('sub-%03d/sub-%03d_corrpowlcmv_gamma',subj,  subj)), 'source', 'pow', 'X', 'tlckpst')%, 'pupild', 'distance');
end
if docorrpow_lcmv_lowfreq
    peakpicking;
    
    datadir = [project_dir 'analysis/'];
    load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_lcmv',  subj,  subj)));
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
    
    load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_source_low',  subj,  subj)));
    load(fullfile(datadir, 'freq/', sprintf('sub-%03d/sub-%03d_freqshort_low', subj, subj)));
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
    
    
    datadir = [project_dir 'analysis/'];
    save(fullfile(datadir, 'corr/', sprintf('sub-%03d/sub-%03d_corrpowlcmv_low', subj, subj)), 'source', 'pow', 'X', 'tlckpst');
end
if docorr_gamma_rt
    load(sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_rt.mat',subj, subj));
    datadir = [project_dir 'analysis/'];
    load(fullfile(datadir, 'source/', sprintf('sub-%03d/sub-%03d_source',  subj,  subj)));
    
    pow      = (abs(F*transpose(freq_shift.fourierspctrm)).^2)*P;
    pow      = log10(pow);
    s = std(pow, [], 2);
    u = mean(pow, 2);
    pow = (pow-repmat(u, [1 size(pow,2)]))./repmat(s, [1 size(pow,2)]);
    
    rho = corr(rt, pow', 'type', 'spearman'); %MvE
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
        erf_osc_datainfo;
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
    erf_osc_datainfo;
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
    
    filename = '/project/3011085.02/analysis/stat_corr_peakpicking_rt2.mat';
    save(filename, 'stat', 'source');
end