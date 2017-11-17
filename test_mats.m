% load some data
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/vol.mat');
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC.mat');
Nsub=1;
headmodel = vol;clear vol;
grad      = dataFIC.grad;

% get the leadfield for an occipital dipole
cfg           = [];
cfg.grad      = ft_convert_units(grad,'cm');
cfg.headmodel = ft_convert_units(headmodel,'cm');
cfg.grid.pos  = [-5 0 6];
cfg.channel   = 'MEG';
leadfield     = ft_prepare_leadfield(cfg);
LF            = leadfield.leadfield{1}(:,3);

% create a 'canonical erf' time course
erf = [zeros(1,100), sin(2.*pi.*(0:99)./100) zeros(1,100)];

% define the relationship between the gamma amplitude and the erf emplitude
ivar  = (1+rand(200,1));
tmp = 1+rand(200,1);
alpha = 0.2;
amp   = alpha.*ivar + (1-alpha).*rand(200,1);
amp2   = alpha.*tmp + (1-alpha).*rand(200,1);

erf_trials2 = amp2*erf(145:155);

% scale the erf time course with the trial specific amplitude
erf_trials = amp*erf;
baseline = 0*amp*ones(1,300);
% for subj=1:Nsub
    keep LF alpha amp cfg dataFIC erf erf_trials grad headmodel ivar leadfield tlckpc tlckpc_bl subj baseline Nsub erf_trials2
    % define some common noise amplitude
    common_noise = (1+rand(200,1));
    amp_noise    = 1e-10.*common_noise;
    amp_noise_bl = 1e-10.*common_noise;
    
    % create single trial data
    trial = zeros(200,151,300);
    trial_bl = trial;
    for k = 1:200
        trial(k,:,:) =  LF*erf_trials(k,:) + randn(151,300).*amp_noise(k);
%                 trial(k,:,145:155) =  LF*erf_trials2(k,:) + randn(151,11).*amp_noise(k);
        trial_bl(k,:,:) = LF*baseline(k,:) + randn(151,300).*amp_noise_bl(k);
    end
%     trial(:,:,145:155)=0;
    
%     shift the occurance of the erf in time
%     shift = randi(21,200,1)-11;
%     for k=1:200
%     trial_shifted(k,:,:) = circshift(trial(k,:,:), shift(k) , 3);
%     end
%     trial=trial_shifted;
    
    % create the design matrix, where the gamma regressor is created by
    % corrupting the ivar with some common noise
    gamma  = ivar + common_noise;
    design = [ones(1,200);gamma'-mean(gamma)];
    
    avg_trial = repmat(mean(trial,3), [1 1 300]);
    avg_trialbl = repmat(mean(trial_bl,3), [1 1 300]);
    % do the regression
    for k = 1:151
        tmp = (design*design')\design*(squeeze(trial(:,k,:))-squeeze(avg_trial(:,k,:)));
        B(k,:) = tmp(2,:);
        tmp = (design*design')\design*(squeeze(trial_bl(:,k,:))-squeeze(avg_trialbl(:,k,:)));
        Bbl(k,:) = tmp(2,:);
    end
    
    
    tlck = [];
    tlck.avg = B;
    tlck.time = (0:299)./100;
    tlck.dimord = 'chan_time';
    tlck.label = leadfield.label;
    tlck.grad = grad;
    
    tlck_bl = [];
    tlck_bl.avg = Bbl;
    tlck_bl.time = (0:299)./100;
    tlck_bl.dimord = 'chan_time';
    tlck_bl.label = leadfield.label;
    tlck_bl.grad = grad;
    
    load ctf151_neighb;
    cfg = [];
    cfg.planarmethod = 'sincos';
    cfg.neighbours = neighbours;
    tlckp = ft_megplanar(cfg,tlck);
    tlckp_bl = ft_megplanar(cfg,tlck_bl);
    
    tlckpc{subj} = ft_combineplanar([],tlckp);
    tlckpc_bl{subj} = ft_combineplanar([],tlckp_bl);
   
    % cfg = [];
    % cfg.operation = 'subtract';
    % cfg.parameter = 'avg';
    % tlckpc_diff = ft_math(cfg, tlckpc, tlckpc_bl);
% end


cfg=[];
cfg.appenddim = 'rpt';
A = ft_appendtimelock(cfg, tlckpc{1:Nsub});
Aavg = ft_timelockanalysis([],A);
B = ft_appendtimelock(cfg, tlckpc_bl{1:Nsub});
Bavg = ft_timelockanalysis([], B);
C = B;
C.trial = repmat(mean(C.trial,3), [1 1 300]);
Cavg = ft_timelockanalysis([], C);

cfg=[];
cfg.parameter = 'avg';
cfg.operation = 'subtract';
d1 = ft_math(cfg, Aavg, Bavg);
d2 = ft_math(cfg, Aavg, Cavg);

cfgp=[];
cfgp.layout = 'CTF151_helmet.mat';
figure;ft_singleplotER(cfgp,tlckpc{1},tlckpc_bl{1});
% figure; ft_singleplotER(cfgp, d1);
% figure; ft_singleplotER(cfgp, d2);
% cfgp.xlim=[1.5 1.5];
% figure; ft_topoplotER(cfgp, d1);

%%

cfg                  = [];
cfg.method           = 'template';
cfg.feedback         = 'no';
neighbours           = ft_prepare_neighbours(cfg, A); % define neighbouring channels

cfg                  = [];
cfg.neighbours       = neighbours;
cfg.parameter        = 'trial';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = 0.05;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.correcttail      = 'prob';
cfg.numrandomization = 1000;


cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg, A, B);
%%
% cfgp = [];
% cfgp.layout = 'CTF151_helmet.mat';
% cfgp.zlim = 'maxabs';
% figure;ft_topoplotER(cfgp,tlckpc_diff);
% figure; ft_singleplotER(cfgp, tlckpc_diff)

% figure;ft_topoplotER(cfgp,tlckpc,tlckpc_bl);
