
if ~exist('doglm_parc', 'var'), doglm_parc = false; end
if ~exist('dosplitpow_source', 'var'), dosplitpow_source = false; end

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
        tmp{k} = load(fullfile([datadir, sprintf('sub-%03d_glm_parcel.mat', subj)]));
        tmp2{k} = load(fullfile([datadir, sprintf('sub-%03d_parcel_blstd.mat', subj)]));
        k=k+1;
    end
    
    for k=1:32
        betas{k} = tmp{k}.betas;
        tlck{k} = tmp{k}.tlck;
        baseline_std{k} = tmp2{k}.baseline_std;
    end
    clear tmp
    
    cfg=[];
    cfg.appenddim = 'rpt';
    tlck_all  = ft_appendtimelock(cfg, tlck{:});
    betas_all  = ft_appendtimelock(cfg, betas{:});
    baseline_std = cat(2, baseline_std{:});
    clear tlck betas
    
    % save average as reference
    cfg=[];
    cfg.avgoverrpt = 'yes';
    tlck_orig = ft_selectdata(cfg, tlck_all);
    
    %% Align parcels over subjects using SVD
    
    for k = 1:length(tlck_all.label);
        erf = squeeze(tlck_all.trial(:,k,:));
        erf = diag(1./std(erf,[],2))*erf; % normalize before SVD
        
        [u,s,v]=svd(erf);
        u2(:,k) = sign(u(:,1));
    end
    u2 = repmat(u2, [1,1, length(tlck_all.time)]);
    
    tlck_all.trial = u2.*tlck_all.trial;
    betas_all.trial = u2.*tlck_all.trial;
    
    cfg=[];
    cfg.avgoverrpt = 'yes';
    tlck_svd1 = ft_selectdata(cfg, tlck_all);
    betas_svd1 = ft_selectdata(cfg, betas_all);

    %% Align ambiguous polarity over parcels using SVD
    if strcmp(erfoi, 'early')
        t1 = 0;
        t2 = 0.2;
    elseif strcmp(erfoi, 'late')
        t1 = 0.2;
        t2 = 0.5;
    end
    t1_svd = nearest(tlck_all.time, t1);
    t2_svd = nearest(tlck_all.time, t2);
    
    erf = squeeze(mean(tlck_all.trial,1));
    [u,s,v]=svd(erf(:,t1_svd:t2_svd),'econ');
    
    tlck_all.trial = permute(repmat(sign(u(:,1)),[1, size(tlck_all.trial,1), length(tlck_all.time)]),[2,1,3]).*tlck_all.trial;
    betas_all.trial = permute(repmat(sign(u(:,1)),[1, size(betas_all.trial,1), length(betas_all.time)]),[2,1,3]).*betas_all.trial;
    
    cfg=[];
    cfg.avgoverrpt = 'yes';
    tlck_svd2 = ft_selectdata(cfg, tlck_all);
    betas_svd2 = ft_selectdata(cfg, betas_all);    


    %% do statistics
    
    ref = betas_all;
    ref.trial = ref.trial*0;
    
    
    cfgs.parameter = 'trial';
    cfgs.clusteralpha = 0.05;
    load parcellation374_neighb
    cfgs.neighbours = neighbours;
    cfgs.latency = [t1 t2];
    cfgs.avgovertime = 'yes';
    stat = ft_timelockstatistics(cfgs, betas_all, ref);
    
    save(fullfile([savedir, sprintf('stat_glm_parcel_%s.mat', erfoi)]), 'stat', 'tlck_orig', 'tlck_svd1', 'tlck_svd2', 'betas_svd2', 'betas_svd1', 't1', 't2');
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