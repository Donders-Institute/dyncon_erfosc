% load some data
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/vol.mat');
load('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataFIC.mat');

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
alpha = 0.2;
amp   = alpha.*ivar + (1-alpha).*rand(200,1); 

% scale the erf time course with the trial specific amplitude
erf_trials = amp*erf;

% define some common noise amplitude
common_noise = (1+rand(200,1));
amp_noise    = 1e-10.*common_noise;
amp_noise_bl = 1e-10.*common_noise;

% define a DC-offset that is correlated with the ivar, with a static
% spatial gradient
grad  = ft_datatype_sens(grad);
[a,b] = match_str(leadfield.label,grad.label);
ypos  = grad.chanpos(b,2);
dc_gradient = -ypos./max(ypos).*5e-10;

% create single trial data 
trial = zeros(200,151,300);
trial_bl = trial;
for k = 1:200
  trial(k,:,:) = LF*erf_trials(k,:) + randn(151,300).*amp_noise(k) + dc_gradient*ones(1,300).*ivar(k);
  trial_bl(k,:,:) = randn(151,300).*amp_noise_bl(k) + dc_gradient*ones(1,300).*ivar(k);
end


% create the design matrix, where the gamma regressor is created by
% corrupting the ivar with some common noise
gamma  = ivar;% + common_noise;
design = [ones(1,200);gamma'-mean(gamma)];

% do the regression
for k = 1:151
  tmp = (design*design')\design*squeeze(trial(:,k,:));
  B(k,:) = tmp(2,:);
  tmp = (design*design')\design*squeeze(trial_bl(:,k,:));
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

tlckpc = ft_combineplanar([],tlckp);
tlckpc_bl = ft_combineplanar([],tlckp_bl);

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
tlckpc_diff = ft_math(cfg, tlckpc, tlckpc_bl);


cfgp = [];
cfgp.layout = 'CTF151_helmet.mat';
figure;ft_topoplotER(cfgp,tlckpc_diff);


figure;ft_topoplotER(cfgp,tlckpc,tlckpc_bl);

% what happens, if we now subtract the baseline average prior to combining
% planar gradients
X = mean(tlckp_bl.avg,2);
%%
cfg           = [];
cfg.parameter = 'avg';
cfg.operation = 'subtract';
cfg.matrix    = repmat(X,[1 300]);
tlckpc2    = ft_combineplanar([], ft_math(cfg, tlckp));
tlckpc2_bl = ft_combineplanar([], ft_math(cfg, tlckp_bl));

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
tlckpc2_diff = ft_math(cfg, tlckpc2, tlckpc2_bl);