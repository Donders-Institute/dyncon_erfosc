%%%%%%%%%%%%%%%%%%%%
% plotting figures %
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% FIGURE 2 %
%%%%%%%%%%%%
%% Single trial ERF (SNR)
subj = 13;

erf_osc_datainfo;
subject = subjects(subj);

[p,f,e]       = fileparts(subject.dataset);
basedir       = strrep(p, 'raw', 'processed');
filename_data = fullfile(basedir, 'cleandata.mat');
load(filename_data);

[data_onset, data_shift] = erfosc_getdata(dataClean, []);
clear dataClean;

% identify channel-of-interest
cfg = [];
cfg.preproc.demean='yes';
cfg.preproc.baselinewindow = [-.1 0];
cfg.vartrllength = 2;
tlck = ft_timelockanalysis(cfg, data_shift);

% use ft_topoplotER, I selected MRP52, based on the 0.07-0.08 time window,
% indx 221

% identify parcel-of-interest
loaddir = '/home/language/jansch/erfosc';
load(fullfile(loaddir, sprintf('sub-%03d_lcmv', subj)), 'source_parc', 'noise');
% I selected R_18_B05_04, indx 360

dat1 = cellrowselect(data_shift.trial, 221);
dat2 = source_parc.F{360}(1,:)*data_shift.trial;

dat1 = cat(1,dat1{:}); dat1 = dat1(:,391:end);
dat2 = cat(1,dat2{:}); dat2 = dat2(:,391:end);

tim = tlck.time(391:end);

figure;imagesc(tim,1:size(dat1,1),ft_preproc_lowpassfilter(dat1,600,40,[],'firws'), [-0.8e-12 0.8e-12]); title('single sensor, filtered'); xlim([0 0.3]); colormap(flipud(brewermap(64, 'RdBu')));
figure;imagesc(tim,1:size(dat2,1),ft_preproc_lowpassfilter(dat2,600,40,[],'firws'),[-3e-12 3e-12]); title('single parcel, filtered'); xlim([0 0.3]); colormap(flipud(brewermap(64, 'RdBu')));



%% ERF topographies

subj = 13;

erf_osc_datainfo;
subject = subjects(subj);

[p,f,e]       = fileparts(subject.dataset);
basedir       = strrep(p, 'raw', 'processed');
filename_data = fullfile(basedir, 'cleandata.mat');
load(filename_data);

[data_onset, data_shift] = erfosc_getdata(dataClean, []);
clear dataClean;

% identify channel-of-interest
cfg = [];
cfg.preproc.demean='yes';
cfg.preproc.baselinewindow = [-.1 0];
cfg.vartrllength = 2;
tlck = ft_timelockanalysis(cfg, data_shift);

% channel
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.xlim=[0.07 0.08];
cfgp.highlight = 'on';
cfgp.highlightchannel = 'MRP52';
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
cfgp.marker = 'off';
cfgp.colorbar='yes';
cfgp.contournum = 0;
cfgp.gridscale = 250;
cfgp.zlim=[-1e-13 1e-13];
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
figure; ft_topoplotER(cfgp, tlck);


% source
loaddir = '/home/language/jansch/erfosc';
load(fullfile(loaddir, sprintf('sub-%03d_lcmv', subj)), 'source_parc', 'noise');
load('atlas_subparc374_8k.mat')
load cortex_inflated_shifted.mat; atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
idx = 1:374;
idx(exclude_label) = [];

  for k = 1:numel(source_parc.label)
    F(k,:) = source_parc.F{k}(1,:);
  end
  dat = F*data_shift.trial;
  
  source=[];
  tim = data_shift.time{1};
  source.avg = zeros(374,1);
  tmp = mean(cat(3,dat{:}),3);
  t1 = nearest(tim, 0.07);
  t2 = nearest(tim, 0.08);
source.avg(idx) = mean(tmp(:,t1:t2),2);
source.avg(exclude_label) = nan;
source.time = 0.075;
source.dimord = 'chan_time';
source.label = atlas.parcellationlabel;
source.brainordinate = atlas;

% show absolute values of ERF (polarities are ambiguous)
source.avg = abs(source.avg);


cfgp = [];
cfgp.method='surface';
cfgp.funparameter='avg';
cfgp.funcolormap = flipud(brewermap(64, 'RdBu'));
cfgp.camlight = 'no';
cfgp.colorbar = 'no';
cfgp.funcolorlim = [-12e-13 12e-13];
cfgp.maskstyle = 'colormix';
cfgp.maskparameter = cfgp.funparameter;
ft_sourceplot(cfgp,source)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;

view([64 16])
view([116 16])
view([-53 2])
view([-120 2])






%%%%%%%%%%%%
% FIGURE 3 %
%%%%%%%%%%%%
%% Time-Frequency plot + Topography (Channel level)

load /project/3011085.02/results/TFR_all.mat

clear all
erf_osc_datainfo;

% subject 26
% TFR
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0 1.5];
cfgp.zlim=[-2.5 2.5];
n=11;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.channel = {'MZO02', 'MRO11', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32'};
cfgp.trials = 25;
ft_singleplotTFR(cfgp, d);

% topography
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0.16 1.6];
cfgp.ylim=[43 56];
cfgp.zlim=[-1.5 1.5];
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MZO02','MRO11', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
n=7;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.trials = 25;
cfgp.colorbar='yes';
cfgp.numcontour = 0;
cfgp.gridscale=250;
ft_topoplotTFR(cfgp, d);


% Subject grand average
% TFR
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0 1.5];
cfgp.zlim=[-1 1];
n=11;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.channel = {'MLO21', 'MLO22', 'MLO31', 'MRO21', 'MRO22', 'MRO31', 'MZO02'};
ft_singleplotTFR(cfgp, d);

% Topography
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0.15 1.05];
cfgp.ylim=[48 64];
cfgp.zlim=[-0.8 0.8];
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MLO21', 'MLO22', 'MLO31', 'MRO21', 'MRO22', 'MRO31', 'MZO02'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
n=9;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.colorbar='yes';
cfgp.numcontour = 0;
cfgp.gridscale = 250;
ft_topoplotTFR(cfgp, d);



%% Gamma-power source plot

erf_osc_datainfo;
sub=26;
mri = ft_read_mri(fullfile([subjects(sub).mridir, '/preproc/mni_resliced.mgz']));
transform = ft_read_mri(fullfile([subjects(sub).mridir, '/preproc/transform_vox2ctf.mat']));
mri.coordsys                = 'ctf';
mri.transform               = transform;

k=1;
for subj=allsubs
tmp = load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel', subj), 'sourceDiff');
sourceDiff{k} = tmp.sourceDiff;
[val, idx] = sort(sourceDiff{k}.avg.pow, 'descend');
idx(isnan(val))=[];
val(isnan(val))=[];
k=k+1;
end


% Grand average
s=sourceDiff{1};
tmp=[];
for k=1:32
tmp(:, k) = sourceDiff{k}.avg.pow;
end
s.avg.pow = mean(tmp,2);

mri_template = ft_read_mri('/project/3011085.02/scripts/fieldtrip/template/anatomy/single_subj_T1_1mm.nii');
load('/project/3011085.02/scripts/fieldtrip/template/sourcemodel/standard_sourcemodel3d6mm.mat');
s.pos = sourcemodel.pos;
cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'avg.pow';
sourceInt  = ft_sourceinterpolate(cfg, s , mri_template);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'pow';
cfg.funcolorlim   = [-1 1];
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap    = 'rampup'; 
cfg.slicerange = [30 35];
cfg.nslices=1;
cfg.renderer = 'painter';
cfg.gridscale = 250;
n=11;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfg.funcolormap=(cmap([2:n n:end-1],:));
ft_sourceplot(cfg, sourceInt);


% subject 26
cfg            = [];
cfg.downsample = 2;
cfg.parameter = 'avg.pow';
sourceInt  = ft_sourceinterpolate(cfg, sourceDiff{25} , mri);

cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.opacitymap    = 'rampup'; 
cfg.funcolorlim   = [-3 3];
cfg.slicerange = [55 60];
cfg.nslices=1;
cfg.renderer = 'painter';
cfg.gridscale = 250;
n=9;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfg.funcolormap=(cmap([1:n n:end],:));
ft_sourceplot(cfg, sourceInt);



%% Gamma-power boxplots
erf_osc_datainfo;

k=1;
for subj=allsubs
tmp{k} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj));
k=k+1;
end


for k=1:32
peakFreq(k) = tmp{k}.peakFreq_gamma;
ratio_broadband(k) = tmp{k}.gamRatio; % mean power increase over channel and frequency
powRatio{k} = tmp{k}.powRatio;
end


for k=1:32
f(k) = find(powRatio{k}.freq==peakFreq(k));
ratio_maxchan(k) = max(powRatio{k}.powspctrm(:,f(k))); % power at max channel at peak frequency
k=k+1;
end

ratio_maxchan_GA = mean(ratio_maxchan);
ratio_maxchan_SD = std(ratio_maxchan);

addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
iosr.statistics.boxPlot(ratio_maxchan'-1, 'showScatter', true, 'scatterMarker', '.')
iosr.statistics.boxPlot(peakFreq', 'showScatter', true, 'scatterMarker', '.')







%%%%%%%%%%%%
% Figure 4 %
%%%%%%%%%%%%
%% Correlation Gamma ERF
load('/project/3011085.02/results/stat_peakpicking3.mat', 'stat','S');
load atlas_subparc374_8k.mat
load cortex_inflated_shifted; atlas.pos=ctx.pos;
stat.brainordinate=atlas;

cfgx = [];
cfgx.method='surface';
cfgx.funparameter='stat';
cfgx.funcolormap = flipud(brewermap(64, 'RdBu'));

cfgx.funcolorlim = [-3 3];
cfgx.maskstyle='colormix';
cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
cfgx.colorbar = 'no';
ft_sourceplot(cfgx,stat)

h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;
view([64 16])
view([116 16])
view([-53 2])
view([-120 2])

% scatter plot correlations per subject. (average correlation over parcels
% that contribute to cluster)
x=find(stat.mask==1);
y = mean(S.rho(:,x),2);
scatter(1:32, sort(y),'.');
addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(y)

%% Correlation ERF-rt

load('/project/3011085.02/results/stat_corr_peakpicking_rt.mat');
load atlas_subparc374_8k.mat
load cortex_inflated_shifted; atlas.pos=ctx.pos;

stat.brainordinate = atlas;
cfgx = [];
cfgx.method='surface';
cfgx.funparameter='stat';
cfgx.funcolormap = flipud(brewermap(64, 'RdBu'));
cfgx.funcolorlim = [-3 3];
cfgx.maskstyle='colormix';
cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
cfgx.colorbar = 'no';
cfgx.opacitylim = 'minzero';
ft_sourceplot(cfgx,stat)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;
view([64 16])
view([130 16])
view([-53 2])
view([-130 2])

R = cat(2,rho{:});
x=find(stat.negclusterslabelmat==1);
y = mean(R(:,x),2);
scatter(1:32, sort(y),'.');
addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(y)


