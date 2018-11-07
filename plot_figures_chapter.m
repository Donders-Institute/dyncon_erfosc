%%%%%%%%%%%%%%%%%%%%
% plotting figures %
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%
% FIGURE 2 %
%%%%%%%%%%%%
%% Gamma Time-Frequency plot + Topography (Channel level)

load /project/3011085.02/results/tfr_all_high_onset.mat

clear all
erf_osc_datainfo;

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

%% Gamma power on 2D sourcemodel (which is used in correlation with ERF)

erf_osc_datainfo;
datadir = '/home/language/jansch/erfosc';
load cortex_inflated_shifted;

k=1;
for sub=allsubs
    if k==1;
        tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source',    sub)), 'source_shift', 'Tval');
        source_shift = tmp{k}.source_shift;
    else
        tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source',    sub)), 'Tval','source_shift');
    end
    Tval(:,k) = tmp{k}.Tval;
    k=k+1;
end

source_shift.pos = ctx.pos;
source_shift.Tval = mean(Tval,2);

cfgx = [];
cfgx.method='surface';
cfgx.funparameter='Tval';
cfgx.funcolorlim = 'maxabs';%[-25 25];
cfgx.maskstyle='colormix';
cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
cfgx.colorbar = 'no';
cfgx.opacitylim = 'zeromax';
n=11;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgx.funcolormap=(cmap([2:n n:end-1],:));
ft_sourceplot(cfgx,source_shift)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;
view([64 16])
view([124 14])
view([-48 5])
view([-120 9])


% include virtual channel positions
datadir = '/home/language/jansch/erfosc';
k=1;
for subj=allsubs;
    tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source',    subj)));
%     load(fullfile(datadir, sprintf('sub-%03d_freqshort',    subj)));
    [m, idx(k)] = max(tmp{k}.Tval);
    k=k+1;
end

grad.chanpos=ctx.pos(idx,:);
grad.elecpos=ctx.pos(idx,:);
for k=1:32
grad.label{k} = sprintf('sub%d', k);
end
ft_plot_sens(grad,'elecshape','sphere', 'elecsize', 8, 'facecolor', 'black')

%% Gamma-power boxplots + power spectra
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
    rat(k,:) = mean(powRatio{k}.powspctrm,1);
end
figure; plot(powRatio{1}.freq, rat, 'gray'); hold on; plot(powRatio{k}.freq, mean(rat,1), 'LineWidth', 4);

for k=1:32
f(k) = find(powRatio{k}.freq==peakFreq(k));
ratio_maxchan(k) = max(powRatio{k}.powspctrm(:,f(k))); % power at max channel at peak frequency
k=k+1;
end

ratio_maxchan_GA = mean(ratio_maxchan);
ratio_maxchan_SD = std(ratio_maxchan);

addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(ratio_maxchan'-1, 'showScatter', true, 'scatterMarker', '.')
figure; iosr.statistics.boxPlot(peakFreq', 'showScatter', true, 'scatterMarker', '.')



%%%%%%%%%%%%
% FIGURE 3 %
%%%%%%%%%%%%
%% low frequency TFR + topo
load /project/3011085.02/results/tfr_all_low_onset.mat

clear all
erf_osc_datainfo;

% Subject grand average
% TFR
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0 1.5];
cfgp.zlim=[-0.3 0.3];
n=7;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.channel = {'MLO23', 'MLO32', 'MLO33', 'MRO23', 'MRO24', 'MRO32', 'MRO33', 'MRO34', 'MRO43', 'MRO44', 'MRO53', 'MRT57'};
ft_singleplotTFR(cfgp, d);

% Topography
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0.5 1.5];
cfgp.ylim=[8 20];
cfgp.zlim=[-0.3 0.3];
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MLO23', 'MLO32', 'MLO33', 'MRO23', 'MRO24', 'MRO32', 'MRO33', 'MRO34', 'MRO43', 'MRO44', 'MRO53', 'MRT57'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
n=7;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.colorbar='yes';
cfgp.numcontour = 0;
cfgp.gridscale = 250;
ft_topoplotTFR(cfgp, d);


%% low frequency power on 2D sourcemodel
erf_osc_datainfo;
datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
load cortex_inflated_shifted;

k=1;
for sub=allsubs
    if k==1;
        tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source_low',    sub)), 'source_shift', 'Tval');
        source_shift = tmp{k}.source_shift;
    else
        tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source_low',    sub)), 'Tval');
    end
    Tval(:,k) = tmp{k}.Tval;
    k=k+1;
end
source_shift.Tval = mean(Tval,2);
source_shift.pos = ctx.pos;

cfgx = [];
cfgx.method='surface';
cfgx.funparameter='Tval';
cfgx.funcolorlim = [-8 8];
cfgx.maskstyle='colormix';
cfgx.maskparameter = cfgx.funparameter;
cfgx.camlight = 'no';
cfgx.colorbar = 'no';
cfgx.opacitylim = 'minzero';
n=9;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgx.funcolormap=(cmap([2:n n:end-1],:));
ft_sourceplot(cfgx,source_shift)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;
view([64 16])
view([119 15])
view([-48 5])
view([-120 9])

% include virtual channel positions
 datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data';
k=1;
for subj=allsubs;
    tmp{k} =     load(fullfile(datadir, sprintf('sub-%03d_source_low',    subj)));
    [m, idx(k)] = min(tmp{k}.Tval);
    k=k+1;
end

grad.chanpos=ctx.pos(idx,:);
grad.elecpos=ctx.pos(idx,:);
for k=1:32
grad.label{k} = sprintf('sub%d', k);
end
ft_plot_sens(grad,'elecshape','sphere', 'elecsize', 8, 'facecolor', 'black')

%% low frequency power boxplots and powerspectra
erf_osc_datainfo;

k=1;
for subj=allsubs
tmp{k} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow_low.mat', subj));
k=k+1;
end

for k=1:32
    spectrum(k,:) = tmp{k}.powratio_occipital;
end
plot(2:2:30, spectrum); hold on; plot(2:2:30, mean(spectrum,1), 'LineWidth', 4);


for k=1:32
peakfreq(k) = tmp{k}.peakfreq;
ratio(k) = tmp{k}.powratio_min; % mean power increase over channel and frequency
end

addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(ratio', 'showScatter', true, 'scatterMarker', '.'); ylim([-0.6 0]);
figure; iosr.statistics.boxPlot(peakfreq', 'showScatter', true, 'scatterMarker', '.')


%%%%%%%%%%%%
% FIGURE 4 %
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
% Figure 5 %
%%%%%%%%%%%%
%% Correlation Gamma ERF
load('/project/3011085.02/results/stat_peakpicking_gamma.mat', 'stat','S');
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
cfgx.opacitylim = 'zeromax';
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

load('/project/3011085.02/results/stat_peakpicking_gamma.mat', 'stat_eye','S_eye');
x2=find(stat_eye.mask==1);
y2 = mean(S_eye.rho(:,x),2);

[val, idx] = sort(y);

figure; scatter(1:32, val,'.'); hold on; scatter(1:32, y2(idx), '.');
addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(y,'showScatter',false,'scatterMarker', '.'); 
figure; iosr.statistics.boxPlot(y2,'showScatter',false,'scatterMarker', '.')
%%%%%%%%%%%%
% Figure 6 %
%%%%%%%%%%%%
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
view([123 16])
view([-53 2])
view([-130 2])

R = source.rho;
x=find(stat.negclusterslabelmat==1);
y = mean(R(:,x),2);
scatter(1:32, sort(y),'.');
addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(y,'showScatter',true,'scatterMarker','.')





%%%%%%%%%%%%%%%%%%%%%%
% Reviewer questions
%%%%%%%%%%%%%%%%%%%%%%
%% reviewer 1, comment 1
% not necessary to plot anything. Absolute coordinates are not saved
% because everything is monitored relative to the start of the experiment.
erf_osc_datainfo;


cfg=[];
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                              'HLC0021','HLC0022','HLC0023', ...
                              'HLC0031','HLC0032','HLC0033'};
k=1;
for subj=allsubs
tmp{k} = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/headmotion.mat', subj));
headpos{k} = ft_selectdata(cfg, tmp{k}.headmotion);
k=k+1;
end
% HLC00n1 X coordinate relative to the dewar (in meters) of the nth head localization coil
% HLC00n2 Y coordinate relative to the dewar (in meters) of the nth head localization coil
% HLC00n3 Z coordinate relative to the dewar (in meters) of the nth head localization coil

% calculate the mean coil position per trial, in x,y,z.
for k=1:32
ntrials = length(headpos{k}.sampleinfo);
for t = 1:ntrials
coil1{k}(:,t) = [mean(headpos{k}.trial{1,t}(1,:)); mean(headpos{k}.trial{1,t}(2,:)); mean(headpos{k}.trial{1,t}(3,:))];
coil2{k}(:,t) = [mean(headpos{k}.trial{1,t}(4,:)); mean(headpos{k}.trial{1,t}(5,:)); mean(headpos{k}.trial{1,t}(6,:))];
coil3{k}(:,t) = [mean(headpos{k}.trial{1,t}(7,:)); mean(headpos{k}.trial{1,t}(8,:)); mean(headpos{k}.trial{1,t}(9,:))];
end

% calculate the headposition and orientation per trial
% pos x y z, orientation x, y, z
cc{k} = circumcenter(coil1{k}, coil2{k}, coil3{k});

coil1avg(k,:) = mean(coil1{k},2)';
coil2avg(k,:) = mean(coil2{k},2)';
coil3avg(k,:) = mean(coil3{k},2)';
end

histogram(coil1avg(:,1), 9, 'BinLimits', [min(coil1avg(:,1)), max(coil1avg(:,1))]);
hold on
histogram(coil1avg(:,2), 9, 'BinLimits', [min(coil1avg(:,2)), max(coil1avg(:,2))]);
histogram(coil1avg(:,3), 9, 'BinLimits', [min(coil1avg(:,3)), max(coil1avg(:,3))]);


%% reviewer 1, comment 3
erf_osc_datainfo;
datadir = '/home/language/jansch/erfosc';
k=1;
for subj=allsubs;
    tmp{k} = load(fullfile(datadir, sprintf('sub-%03d_source',    subj)));
%     load(fullfile(datadir, sprintf('sub-%03d_freqshort',    subj)));
    [m, idx(k)] = max(tmp{k}.Tval);
    k=k+1;
end
%   load cortex_inflated_shifted.mat 
% % load cortex_inflated.mat
% % load atlas_subparc374_8k.mat
% % ctx.pos=atlas.pos;
% ctx.vchan = zeros(size(ctx.brainstructure));
% ctx.vchan(idx) = true;
% 
% cfg=[];
% cfg.method = 'surface';
% cfg.funparameter = 'vchan';
% cfg.maskparameter = 'vchan';
% cfg.maskstyle = 'colormix';
% cfg.camlight = 'yes';
% cfg.funcolormap = brewermap(4, 'Reds');
% ft_sourceplot(cfg, ctx); material dull;
% view([106 18]); lighting gouraud; h = camlight(108,-29, 'infinite');% left
% ft_sourceplot(cfg, ctx); material dull;
% view([68 -12]); lighting gouraud;  h = light('position', [0 -1 0]); % right
% ft_sourceplot(cfg, ctx); material dull;
% view([-65 3]); lighting gouraud; h = camlight(-52,-20, 'infinite'); % right medial
% ft_sourceplot(cfg, ctx); material dull;
% view([-107 -20]); lighting gouraud; h = camlight(-54,4, 'infinite'); % left medial

% or
load cortex_inflated.mat 
load atlas_subparc374_8k.mat
ctx.pos=atlas.pos;
% figure; ft_plot_mesh(ctx); hold on
ctx.dum = zeros(size(ctx.brainstructure));
cfg=[];
cfg.funparameter = 'dum';
cfg.method = 'surface';
cfg.maskparameter = 'dum';
cfg.maskstyle = 'colormix';
cfg.funcolorlim = [0 1];
ft_sourceplot(cfg, ctx); material dull; view(0,9); camlight


grad.chanpos=ctx.pos(idx,:);
grad.elecpos=ctx.pos(idx,:);
for k=1:32
grad.label{k} = sprintf('sub%d', k);
end
ft_plot_sens(grad,'coilshape','sphere')

%% reviewer 2, comment 1
erf_osc_datainfo;
k=1;
for subj=allsubs
trl{k} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/trialinfo.mat', subj));
trl{k} = [trl{k}.trialinfo(:,4), (trl{k}.trialinfo(:,9)-trl{k}.trialinfo(:,8))/1200];
rt{k} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj));
rt{k} = rt{k}.rt;
k=k+1;
end

% find trial number of each rt.
for k=1:32
    tmp = trl{k};
    for j=1:numel(rt{k});
    tmpidx = find(tmp(:,2)==rt{k}(j),1);
    trlidx{k}(j) = tmp(tmpidx,1);
    tmp(1:tmpidx,:) = [];
    end
end

% correlation between trial number and reaction time
for k=1:32
    r(k) = corr(rt{k}, trlidx{k}', 'type','spearman');
end
    u = mean(r); % negative correlation, so people actually become faster throughout the experiment.
    [H,P,~,stats] = ttest(r);
    
%     % how about when we look for every block?
%     blockSize = 40;
%     for k=1:32
%         l=1;
%         for b=1:40:numel(rt{k});
%             idx = find(trlidx{k}>=b & trlidx{k}<b+blockSize);
%             tmpr(1,l) = corr(rt{k}(idx), trlidx{k}(idx)');
%             l=l+1;
%         end
%         rblock(k,1) = mean(tmpr); % average correlation over blocks
%     end
%     ublock = mean(rblock); % positive correlation, people get worse over the block
%     [Hb,Pb,~,statsb] = ttest(rblock);
    
%% reviewer 2, comment 2
erf_osc_datainfo;

% get power at maximum channel and peak frequency.
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
end

ratio_maxchan_GA = mean(ratio_maxchan);

% get reaction times
k=1;
for subj=allsubs
    load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj));
    RT{k} = rt;
    clear rt;
    k=k+1;
end

for k=1:32
    rt(k,1)=mean(RT{k});
end

% correlations with age.
load('/project/3011085.02/age.mat');
age(badsubjects)=[]; % get rid of the bad subject

[r pval] = corr(age, rt)
[r pval] = corr(age, ratio_maxchan')
[r pval] = corr(age, peakFreq')
        

%% reviewer 2, comment 3
erf_osc_datainfo;
k=1;
for subj=allsubs
trl{k} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/trialinfo.mat', subj));
k=k+1;
end

for k=1:32
    ntrials(k,1) = size(trl{k}.trialinfo,1);
    ntrials_noresp(k,1) = sum(trl{k}.trialinfo(:,9)==0 & trl{k}.trialinfo(:,8)~=0);
    ntrials_lateresp(k,1) = sum(trl{k}.trialinfo(:,8)~=0 & (trl{k}.trialinfo(:,9)-trl{k}.trialinfo(:,8))/1200>=0.7);
end

deleted_based_on_resp = (ntrials_noresp./ntrials)*100; % percentage
u = mean(deleted_based_on_resp);




load atlas_subparc374_8k.mat
cd /project/3011085.02/scripts/erfosc/analysis_JM_data
d = dir('*corrpowlcmv_gamma_alltrials.mat');
for k = 1:32
load(d(k).name,'source');
if k==1
S=source;
S.rho = source.rho;
else
S.rho(:,k)=source.rho;
end
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
p = stat.posclusters(1);

u = mean(mean(S.rho(:,stat.mask),2));

%% Reviewer 2, comment 6
erf_osc_datainfo;
k=1;
for subj=allsubs
    load(sprintf('/project/3011085.02/scripts/erfosc/analysis_JM_data/sub-%03d_corr_3Dgamma_rt.mat', subj));
    r(:,k) = rho;
    k=k+1;
    clear rho
end

rho = nanmean(r,2);
load cortex_inflated_shifted.mat 
% load atlas_subparc374_8k.mat
% ctx.pos=atlas.pos;
ctx.rho = rho;
cfg=[];
cfg.funparameter = 'rho';
cfg.funcolormap = flipud(brewermap(64, 'RdBu'));
cfg.funcolorlim = [-0.07 0.07];
cfg.method = 'surface';
cfg.maskparameter = cfg.funparameter;
cfg.maskstyle = 'colormix';
cfg.camlight = 'no';
ft_sourceplot(cfg, ctx); material dull;

h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
view([64 16])
view([116 16])
view([-53 2])
view([-120 2])



