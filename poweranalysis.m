%% poweranalysis

%% load data
erf_osc_datainfo;
data1 = load(sprintf('/project/3011085.02/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', 4), 'dataClean');
load(fullfile([pilotsubjects(subj).segmentedmri, '.mat']));

data1 = data1.dataClean;
fs = 1200;

% select only shift trials, with valid response
idxM        = find(data1.trialinfo(:,5)>0 & data1.trialinfo(:,6)>0);
cfg         = [];
cfg.trials  = idxM;
cfg.channel = 'MEG';
data1        = ft_selectdata(cfg, data1);
%%
for ii=1:20
    cfg        = [];
    cfg.trials = randsample(length(data1.trial),0.5*length(data1.trial));
    data = ft_selectdata(cfg, data1);
    
    % redifine zero point to shift onset
    cfg         = [];
    cfg.offset  = -(data.trialinfo(:,5)-data.trialinfo(:,4));
    dataShift   = ft_redefinetrial(cfg, data);
    
    cfg         = [];
    cfg.latency = [-0.05 0.75];
    dataShift2   = ft_selectdata(cfg, dataShift);
    
    
    % aseo
    s.X = 1; % 1 source signal estimate?
    nComponent = 1;
    
    % run a dss decomposition
    params          = [];
    params.time     = dataShift2.time;
    params.demean   = 'prezero';
    params.pre      = 0.05*fs; % 50 ms pre stim reversal as baseline
    params.pst      = 0.75*fs; % 750 ms after stim reversal
    
    cfg                     = [];
    cfg.method              = 'dss';
    cfg.dss.denf.function   = 'denoise_avg2';
    cfg.dss.denf.params     = params;
    cfg.dss.wdim            = 100;
    cfg.numcomponent        = nComponent;
    cfg.cellmode            = 'yes';
    
    % cfg for plotting
    cfgp            = [];
    cfgp.layout     = 'CTF275_helmet.mat';
    cfgp.component  = 1:nComponent;
    comp            = ft_componentanalysis(cfg, dataShift2);
    figure;ft_topoplotIC(cfgp, comp);drawnow;
    
    % get the 'average'
    [~,~,avgcomp] = denoise_avg2(params,comp.trial,s);
    
    
    h=figure; plot(avgcomp(1,:),'r');
    saveas(h,sprintf('FIG%d.fig',ii));
    waveformInitSet = input('waveformInitSet')

%     waveformInitSet{ii} = [123 161;162 226;227 522 ]';
    waveformInitSet = waveformInitSet(:);
    jitter = input('jitter');
    ASEOiteration = 1;
    [q1, q2] = doASEO(comp, 'waveformInitSet', waveformInitSet, 'jitter', jitter, 'numiteration', ASEOiteration);
    
    % gamma peak
    cfg          = [];
    cfg.channel  = {'MZO', 'MZP', 'MLO', 'MLP', 'MRO', 'MRP'};
    cfg.latency  = [-1+1/fs 0];
    dataBaseline = ft_selectdata(cfg, data);
    cfg.latency  = [0.4+1/fs 1.75]; % take active time window after first erfs
    dataActive   = ft_selectdata(cfg, data);
    
    cfg            = [];
    cfg.foi        = 30:1:100;
    cfg.method     = 'mtmfft';
    cfg.output     = 'pow';
    cfg.tapsmofrq  = 1;
    cfg.taper      = 'dpss';
    cfg.keeptrials = 'no';
    cfg.pad        = 2;
    powActive      = ft_freqanalysis(cfg, dataActive);
    powBaseline    = ft_freqanalysis(cfg, dataBaseline);
    
    cfg           = [];
    cfg.operation = '(x1-x2)./x2';
    cfg.parameter = 'powspctrm';
    powDiff       = ft_math(cfg, powActive, powBaseline);
    
    % average over channels, take the freq with max gamma pow diff
    gammaAvg      = mean(powDiff.powspctrm,1);
    [maxP maxIdx] = max(gammaAvg);
    peakFreq(ii)      = powDiff.freq(maxIdx);
    
    %  gamma virtualchan
    cfg         = [];
    cfg.latency = [-1+1/fs 0];
    dataPre     = ft_selectdata(cfg, data);
    % take second preceding shift (NOTE: is it confounding if this includes
    % grating onset, which has higher gamma peak freq?)
    dataPost    = ft_selectdata(cfg, dataShift);
    dataAll     = ft_appenddata([], dataPost, dataPre);
    
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = 6;
    cfg.foilim    = [peakFreq(ii) peakFreq(ii)];
    freqAll       = ft_freqanalysis(cfg, dataAll);
    
    % calculate power pre and post stimulus
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'powandcsd';
    cfg.tapsmofrq = 6;
    cfg.foilim    = [peakFreq(ii) peakFreq(ii)];
    freqPre       = ft_freqanalysis(cfg, dataPre);
    freqPost      = ft_freqanalysis(cfg, dataPost);
    
    % Source analysis
    %{
    % constructs a volume conduction model from the geometry of the head.
    cfg         = [];
    cfg.method  = 'singleshell';
    headmodel   = ft_prepare_headmodel(cfg, segmentedmri);
    headmodel   = ft_convert_units(headmodel, 'cm');
    
    % constructs a source model, for example a 3-D grid or a cortical sheet.
    cfg                 = [];
    cfg.grid.resolution = 0.5;
    cfg.grid.unit       = 'cm';
    cfg.headmodel       = headmodel;
    cfg.grad            = data.grad;
    sourcemodel         = ft_prepare_sourcemodel(cfg);
    
    % computes the forward model for many dipole locations on a regular 2D or
    % 3D grid and stores it for efficient inverse modelling
    cfg                 = [];
    cfg.grad            = data.grad;
    cfg.headmodel       = headmodel;
    cfg.reducerank      = 2;
    cfg.grid.resolution = 0.5;   % use a 3-D grid with a 1 cm resolution
    cfg.grid.unit       = 'cm';
    sourcemodel_lf      = ft_prepare_leadfield(cfg);
    %}
    % calculate common filter
    cfg                   = [];
    cfg.grad              = freqAll.grad;
    cfg.method            = 'dics';
    cfg.frequency         = peakFreq(ii);
    cfg.grid              = sourcemodel_lf;
    cfg.headmodel         = headmodel;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.lambda       = '5%';
    cfg.dics.keepfilter   = 'yes';
    cfg.dics.realfilter   = 'yes';
    sourceAll             = ft_sourceanalysis(cfg, freqAll);
    
    cfg.grid.filter       = sourceAll.avg.filter;
    sourcePre             = ft_sourceanalysis(cfg, freqPre);
    sourcePost            = ft_sourceanalysis(cfg, freqPost);
    sourceDiff            = sourcePost;
    sourceDiff.avg.pow    = (sourcePost.avg.pow - sourcePre.avg.pow) ./ sourcePre.avg.pow;
    
    % Create virtual channel of maximal gamma power
    % the following position contains the max gamma power difference
    [maxval, maxpowindx] = max(sourceDiff.avg.pow);
    sourceDiff.pos(maxpowindx, :);
    % we will create a virtual channel based on this location. In order to do
    % this, we have to do timelockanalysis and use an LCMV beamformer, because
    % this will pass the activity at the location of interest with unit gain,
    % while optimally reducing activity from other locations. This filter then
    % can be applied to the original MEG data
    
    cfg                   = [];
    cfg.covariance        = 'yes';
    cfg.vartrllength      = 2;
    cfg.covariancewindow  = 'all';
    tlock                 = ft_timelockanalysis(cfg, dataAll);
    
    cfg                 = [];
    cfg.method          = 'lcmv';
    cfg.headmodel       = headmodel;
    cfg.grid.pos        = sourcemodel.pos(maxpowindx, :);
    cfg.grid.unit       = sourcemodel.unit;
    cfg.lcmv.keepfilter = 'yes';
    cfg.projectmom      = 'yes';
    cfg.fixedori        = 'yes';
    source_idx          = ft_sourceanalysis(cfg, tlock);
    
    beamformerGamPow = source_idx.avg.filter;
    
    gamPowData              = data;
    gamPowData.label        = {'gam_pow'};
    gamPowData.trial        = [];
    for i=1:length(dataShift.trial)
        gamPowData.trial{i} = beamformerGamPow{1} * data.trial{i};
    end
    
    % gamma pow
    cfg                = [];
    cfg.offset         = -(gamPowData.trialinfo(:,5)-gamPowData.trialinfo(:,4));
    gamPowDataShift    = ft_redefinetrial(cfg, gamPowData);
    
    cfg          = [];
    cfg.latency  = [-0.5+1/fs 0];
    dataBl      = ft_selectdata(cfg, gamPowData);
    dataPre     = ft_selectdata(cfg, gamPowDataShift);
    
    peakFreq(ii) = 2*round(peakFreq(ii)/2);
    smoothing = 6;
    % gamma power
    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
    cfg.tapsmofrq   = smoothing;
    cfg.foilim      = [(peakFreq(ii) - 6*smoothing) (peakFreq(ii) + 6*smoothing)];
    cfg.keeptrials  = 'no'; % average baseline over trials
    gamPowPre       = ft_freqanalysis(cfg, dataBl);
    cfg.keeptrials  = 'yes';
    gamPowPost      = ft_freqanalysis(cfg, dataPre);
    
    gamPow = gamPowPost;
    gamPowPre.powspctrm = repmat(gamPowPre.powspctrm, [size(gamPow.powspctrm,1), 1]);
    
    gamPow.powspctrm = (squeeze(gamPow.powspctrm) - gamPowPre.powspctrm)./gamPowPre.powspctrm;
    
    
    % correlation
    
    peakIdx = find(gamPow.freq==peakFreq(ii));
    [rAmp{ii} pAmp{ii}] = corr(gamPow.powspctrm(:,peakIdx), q1.params.amplitude, 'type', 'spearman');
    [rLat{ii} pLat{ii}] = corr(gamPow.powspctrm(:,peakIdx), q1.params.latency, 'type', 'spearman');
end

for k=1:20
    x(k) = rAmp{k}(1);
end
sd = sqrt(20)*std(x); % corrected for bootstrapping
avg = mean(x);
ntest = 4; % t-test for correlation of gamma power with 2 peaks (amplitude and latency)
nsample = sampsizepwr('t', [0 sd], avg, 0.9, [], 0.05/ntest); 