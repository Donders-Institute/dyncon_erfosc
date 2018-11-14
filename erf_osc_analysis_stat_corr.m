function erf_osc_analysis_stat_corr(freqRange, zeropoint, erfoi, compareQuartile, correlation)

if nargin<1 || isempty(freqRange)
    freqRange = 'high';
end
if nargin<2 || isempty(zeropoint)
    zeropoint = 'reversal';
end
if nargin<3 || isempty(erfoi)
    erfoi = 'reversal';
end
if nargin<4 || isempty(compareQuartile)
    compareQuartile = true;
end
if nargin<5 || isempty(correlation)
    correlation = 'gamma_erf_virtualchan';
end

% initiate diary
ft_diary('on')

%% load data
erf_osc_datainfo;

for subj=allsubs
    if ~compareQuartile
        tmp{subj} = load(sprintf('/project/3011085.02/analysis/corr/sub-%03d/sub-%03d_corr_amp_tfr.mat', subj, subj));
        z_act{subj} = tmp{subj}.z_act;
        z_bl{subj} = tmp{subj}.z_bl;
    else
        if strcmp(correlation, 'amp_tfr')
            tmp{subj} = load(sprintf('/project/3011085.02/analysis/corr/sub-%03d/sub-%03d_corr_amp_tfr_quartile.mat', subj, subj));
            q1{subj} = tmp{subj}.tfa_q1;
            q4{subj} = tmp{subj}.tfa_q4;
        elseif strcmp(correlation, 'gamma_erf_virtualchan')
            tmp{subj} = load(sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_erf_virtualchan_%s.mat', subj, subj, erfoi));
            q1{subj} = tmp{subj}.tl_q1;
            q4{subj} = tmp{subj}.tl_q4;
            tlck{subj} = ft_timelockanalysis([], tmp{subj}.tlck);
        end
            
    end
end
clear tmp


%% Baseline correct
% average baseline over time and repeat over time
if ~compareQuartile
    for subj = allsubs
        cfg=[];
        cfg.avgovertime = 'yes';
        z_bl{subj} = ft_selectdata(cfg, z_bl{subj});
        z_bl{subj}.powspctrm = repmat(z_bl{subj}.powspctrm, [1,1, length(z_act{1}.time)]);
        z_bl{subj}.time = z_act{subj}.time;
    end
    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract'; %relative change
    for subj=allsubs
        diffZ{subj} = ft_math(cfg, z_act{subj}, z_bl{subj});
    end
    
    % get grand average
    cfg               = [];
    cfg.appenddim     = 'rpt';
    diffZ_GA = ft_appendfreq(cfg, diffZ{allsubs});
    z_act_GA = ft_appendfreq(cfg, z_act{allsubs});
    z_bl_GA = ft_appendfreq(cfg, z_bl{allsubs});
    
else
    cfg=[];
    cfg.appenddim = 'rpt';
    q1_GA = ft_appendtimelock(cfg, q1{allsubs});
    q4_GA = ft_appendtimelock(cfg, q4{allsubs});
    tlck_GA = ft_appendtimelock(cfg, tlck{allsubs});
    
    if strcmp(correlation, 'gamma_erf_virtualchan')
        % make sure there is no polarity ambiguity by doing SVD on the
        % group data and multiply the subject data by the respective sign
        % of the first SV component.
        erf = squeeze(tlck_GA.trial);
        erf = diag(1./std(erf,[],2))*erf; % normalize before SVD
        
        [u,s,v]=svd(erf);
        
        u2 = repmat(sign(u(:,1)), [1, 1, length(tlck_GA.time)]);
        tlck_GA.trial = u2.*tlck_GA.trial;
        q1_GA.trial = u2.*q1_GA.trial;
        q4_GA.trial = u2.*q4_GA.trial;
    end
    
    
    cfg=[];
    if strcmp(correlation, 'amp_tfr')
        cfg.parameter = 'powspctrm';
    elseif strcmp(correlation, 'gamma_erf_virtualchan')
        cfg.parameter = 'trial';
    end
    cfg.operation = 'subtract';
    diff_GA = ft_math(cfg, q1_GA, q4_GA);
end



%% Do statistics

Nsub = length(allsubs);

cfg                  = [];
cfg.method           = 'template';
cfg.feedback         = 'no';
if ~compareQuartile
    neighbours           = ft_prepare_neighbours(cfg, z_act_GA); % define neighbouring channels
else
    neighbours           = ft_prepare_neighbours(cfg, q1_GA); % define neighbouring channels
end

cfg                  = [];
% cfg.neighbours       = neighbours;
cfg.channel          = 'gam_pow';
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.numrandomization = 10000;


cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

if ~compareQuartile
    stat = ft_freqstatistics(cfg, z_act_GA, z_bl_GA);
else
    stat = ft_freqstatistics(cfg, q1_GA, q4_GA);
    stat.cfg = rmfield(stat.cfg, 'previous');
end



%% save
if ~compareQuartile
    filename = '/project/3011085.02/analysis/stat_corr_amp_tfr';
    save(fullfile([filename '.mat']), 'stat', 'diffZ_GA', 'z_bl_GA', 'z_act_GA', '-v7.3');
else
    filename = '/project/3011085.02/analysis/stat_corr_amp_tfr_quartile';
    save(fullfile([filename '.mat']), 'stat', '-v7.3');
end
ft_diary('off')



