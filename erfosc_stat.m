
if ~exist('doglm_parc', 'var'), doglm_parc = false; end
if ~exist('dosplitpow_source', 'var'), dosplitpow_source = false; end
if ~exist('doresplocked', 'var'), doresplocked = false; end

ft_diary('on')
datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
savedir = '/project/3011085.02/results/';

cfgs=[];
cfgs.statistic        = 'ft_statfun_depsamplesT';
cfgs.design(1,1:64) = [ones(1,32) 2*ones(1,32)];
cfgs.design(2,1:64) = [1:32 1:32];
cfgs.ivar=1;
cfgs.uvar=2;
cfgs.alpha = 0.05;
cfgs.method = 'montecarlo';
cfgs.clusterstatistic = 'maxsum';
cfgs.correctm = 'cluster';
cfgs.numrandomization = 10000;
cfgs.correcttail = 'prob';

if doglm_parc
    erf_osc_datainfo;
    
    k = 1;
    for subj = allsubs
        if doresplocked
            tmp{k} = load(fullfile([datadir, sprintf('sub-%03d_glm_parcelresp.mat', subj)]));
        else
            tmp{k} = load(fullfile([datadir, sprintf('sub-%03d_glm_parcel.mat', subj)]));
        end
    tmp2{k} = load(fullfile([datadir, sprintf('sub-%03d_parcel_blstd.mat', subj)])); %baseline correction not needed; was already done before GLM.
    tmp3{k} = tmp2{k}.baseline_std;
    k=k+1;
    end
    
    for k=1:32
        betas{k} = tmp{k}.betas;
        tlck{k} = tmp{k}.tlck;
    end
    
    
    cfg=[];
    cfg.appenddim = 'rpt';
    tlck_all  = ft_appendtimelock(cfg, tlck{:});
    betas_all  = ft_appendtimelock(cfg, betas{:});
    
    if ~doresplocked
        baselinestd = permute(cat(2, tmp3{:}),[2,1]); % standard deviation over concatenated single trial baselines BEFORE any SVD
        tb1 = nearest(tlck_all.time, -.1);
        tb2 = nearest(tlck_all.time, 0);
        mu_erf = mean(tlck_all.trial(:,:,tb1:tb2),3);
        mu_betas_all = mean(betas_all.trial(:,:,tb1:tb2),3);    
    end
%     clear tlck betas tmp*
    
    % realign the betas in the direction of the effect. If the ERF is
    % positive and the beta as well, the effect is positive. The same goes
    % for a negative ERF and negative beta. If the sign of the ERF and of
    % the betas is opposite, the effect is negative (stronger gamma power
    % leads so weaker ERF).
    cfg=[];
    cfg.parameter = 'trial';
    cfg.operation = 'sign(x1).*x2';
    betas_align = ft_math(cfg, tlck_all, betas_all);
    
    if ~doresplocked
        mu_betas_align = mean(betas_align.trial(:,:,tb1:tb2),3);

        tlck_norm = tlck_all;
        tlck_norm.trial = (tlck_norm.trial.*repmat(mu_erf, [1,1, length(tlck_norm.time)]))./repmat(baselinestd, [1,1,length(tlck_norm.time)]);
        betas_all_norm = betas_all;
        betas_all_norm.trial = (betas_all_norm.trial.*repmat(mu_betas_all, [1,1, length(betas_all_norm.time)]))./repmat(baselinestd, [1,1,length(betas_all_norm.time)]);
        betas_align_norm = betas_align;
        betas_align_norm.trial = (betas_align_norm.trial.*repmat(mu_betas_align, [1,1, length(betas_align_norm.time)]))./repmat(baselinestd, [1,1,length(betas_align_norm.time)]);
    end
    
    %% do statistics
    
    ref = betas_align;
    ref.trial = ref.trial*0;
    
    
    cfgs.parameter = 'trial';
    cfgs.clusteralpha = 0.05;
    load parcellation374_neighb
    cfgs.neighbours = neighbours;
    if doresplocked
        cfgs.latency = [-0.5 0.4];
    else
        cfgs.latency = [0 0.5];
    end
    stat = ft_timelockstatistics(cfgs, betas_align, ref);
    
    if doresplocked
        save(fullfile([savedir, 'stat_glm_parcelresp.mat']), 'stat', 'tlck_all', 'betas_all', 'betas_align');
    else
        save(fullfile([savedir, 'stat_glm_parcel.mat']), 'stat', 'tlck_all', 'betas_all', 'betas_align', 'tlck_norm', 'betas_all_norm', 'betas_align_norm', 'mu_erf', 'mu_betas_all', 'mu_betas_align');
    end
    ft_diary('off')
    
elseif dosplitpow_source
    erf_osc_datainfo;
    for subj=allsubs
        tmp1 = load(fullfile([datadir, sprintf('sub-%03d_splitpow_source_low.mat', subj)]), 'tfr_1', 'tfr_2');
        tfrl_1{subj} = tmp1.tfr_1;
        tfrl_2{subj} = tmp1.tfr_2;
        tmp2 = load(fullfile([datadir, sprintf('/16 Hz smoothing/sub-%03d_splitpow_source_high.mat', subj)]), 'tfr_1', 'tfr_2');
        tfrh_1{subj} = tmp2.tfr_1;
        tfrh_2{subj} = tmp2.tfr_2;
    end
    
    cfg=[];
    cfg.appenddim = 'rpt';
    TFRl_1 = ft_appendfreq(cfg, tfrl_1{allsubs});
    TFRl_2 = ft_appendfreq(cfg, tfrl_2{allsubs});
    TFRh_1 = ft_appendfreq(cfg, tfrh_1{allsubs});
    TFRh_2 = ft_appendfreq(cfg, tfrh_2{allsubs});
    
    cfgs.parameter = 'powspctrm';
    cfgs.latency = [-0.4 0];
    cfgs.clusteralpha = 0.01;
    load parcellation374_neighb
    cfgs.neighbours = neighbours;
    
    stat_l = ft_freqstatistics(cfgs, TFRl_2, TFRl_1);
    stat_h = ft_freqstatistics(cfgs, TFRh_2, TFRh_1);
    
    save(fullfile([savedir, 'stat_splitpow_source']), 'stat_l', 'stat_h');
end