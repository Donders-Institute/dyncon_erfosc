%% plot_figures

%% rt_all
erf_osc_datainfo
k=1;
for subj=allsubs
tmp{k} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj));
k=k+1;
end
for k=1:32
rt{k} = tmp{k}.rt;
end
for k=1:32
med(k) = median(rt{k});
mn(k) = mean(rt{k});
end
histogram(rt)
vline(mean(med), 'g')
vline(mean(mn), 'r')
title('average mean reaction time 371ms; median 356ms')


%% ERF_pospole_GA

%{
load('/project/3011085.02/results/tlck_reversal_all.mat')

% ERF_all
cfg=[];
cfg.latency = [0 0.5];
cfg.channel = {'MRO22', 'MRO23', 'MRO32', 'MRO33', 'MRO34', 'MRO44', 'MRT47'};
tlck_all_pospole = ft_selectdata(cfg, tlck_all);
plot(tlck_all_pospole.time, squeeze(mean(tlck_all_pospole.trial,2)))
hold on
plot(tlck_all_pospole.time, mean(squeeze(mean(tlck_all_pospole.trial,2))), 'LineWidth', 4)
xlim([0 0.3]);ylim([-2.5e-13 2.5e-13])
vline(0.097);

% ERF_topo_GA
cfgp=[];
cfgp.zlim=[-8e-14 8e-14];
cfgp.xlim= [0.097 0.097];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MRO22', 'MRO23', 'MRO32', 'MRO33', 'MRO34', 'MRO44', 'MRT47'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
n=9;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.colorbar='yes';
cfgp.numcontour = 0;
cfgp.gridscale = 250;
ft_topoplotER(cfgp, ft_timelockanalysis([], tlck_all));
%}


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
% n=11;
% cmap = flipud(brewermap(2*n-1,'RdBu'));
% cfgp.colormap=(cmap([2:n n:end-1],:));
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

source2=source;
source2.avg = abs(source2.avg);


cfgp = [];
cfgp.method='surface';
cfgp.funparameter='avg';
cfgp.funcolormap = flipud(brewermap(64, 'RdBu'));
% cfgp.funcolorlim = [-3 3];
% cfgp.maskstyle='colormix';
% cfgp.maskparameter = cfgx.funparameter;
cfgp.camlight = 'no';
cfgp.colorbar = 'no';
cfgp.funcolorlim = [-12e-13 12e-13];
cfgp.maskstyle = 'colormix';
cfgp.maskparameter = cfgp.funparameter;
ft_sourceplot(cfgp,source2)
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;
view([64 16])
view([116 16])
view([-53 2])
view([-120 2])

%{
% ERF_topo_26
cfgp=[];
cfgp.zlim=[-2e-13 2e-13];
cfgp.xlim= [0.075 0.075];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MLO24', 'MLO33', 'MLO34', 'MLO44','MLO43', 'MLT37', 'MLT47'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
cfgp.trials = 25;
n=9;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.colorbar='yes';
cfgp.numcontour = 0;
cfgp.gridscale = 250;
ft_topoplotER(cfgp, tlck_all);
cfgp=[];
cfgp.channel = {'MLO24', 'MLO33', 'MLO34', 'MLO44','MLO43', 'MLT37', 'MLT47'};
cfgp.trials=25;
cfgp.xlim=[0 0.5];
ft_singleplotER(cfgp, tlck_all);
vline(0.075); xlim([0 0.3])%}
%}

%% Imagesc
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

% figure;imagesc(tim,1:size(dat1,1),dat1, [-10e-13 10e-13]); title('single sensor, unfiltered'); colormap(flipud(brewermap(64, 'RdBu')));
% figure;imagesc(tim,1:size(dat2,1),dat2, [-4e-12 4e-12]); title('single parcel, unfiltered');colormap(flipud(brewermap(64, 'RdBu')));
%figure;plot(tim,dat1,'b');hold on;plot(tim,mean(dat1,1),'r','linewidth',2)
%figure;plot(tim,dat2,'b');hold on;plot(tim,mean(dat2,1),'r','linewidth',2)

figure;imagesc(tim,1:size(dat1,1),ft_preproc_lowpassfilter(dat1,600,40,[],'firws'), [-0.8e-12 0.8e-12]); title('single sensor, filtered'); xlim([0 0.3]); colormap(flipud(brewermap(64, 'RdBu')));
figure;imagesc(tim,1:size(dat2,1),ft_preproc_lowpassfilter(dat2,600,40,[],'firws'),[-3e-12 3e-12]); title('single parcel, filtered'); xlim([0 0.3]); colormap(flipud(brewermap(64, 'RdBu')));
%{
subj=2;
% L_17_B05_03
datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
tmp = load(fullfile([datadir, sprintf('sub-%03d_erfparc.mat', subj)]));
tmp = tmp.parceldata_shift;

cfg=[];
cfg.channel = 'L_17_B05*';
cfg.keeptrials='yes';
cfg.vartrllength = 2;
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 30;
cfg.preproc.lpfilttype = 'firws';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
tlck = ft_timelockanalysis(cfg, tmp);
c = find(strcmp('L_17_B05_03', tlck.label));
figure;imagesc(tlck.time, 1:numel(tlck.trialinfo), squeeze((tlck.trial(:,c,:))))%, [-0.6e-11 0.6e-11]);
xlim([0 0.3])
% colormap(flipud(brewermap(64,'RdBu')));
% n=7;
% cmap = flipud(brewermap(2*n-1,'RdBu'));
% colormap(cmap([2:n n:end-1],:));

%MLO34
d = load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat',subj));
cfg=[];
cfg.keeptrials='yes';
cfg.vartrllength = 2;
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 30;
cfg.preproc.lpfilttype = 'firws';
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
chan = ft_timelockanalysis(cfg, d.dataClean);
x=tlck.trialinfo;
for k=1:numel(x)
y(k) = find(chan.trialinfo(:,1)==x(k));
end

cfgp.layout = 'CTF275_helmet.mat';
figure; ft_singleplotER(cfgp, chan);

cfg=[];
cfg.trials=y;
cfg.channel = 'MLO22'%'MLO34';
chan=ft_selectdata(cfg, chan);

figure; imagesc(chan.time, 1:numel(x), squeeze(chan.trial))%, [-5e-13 5e-13])
xlim([0 0.3])
% colormap(flipud(brewermap(64,'RdBu')));
% n=11;
% cmap = flipud(brewermap(2*n-1,'RdBu'));
% colormap(cmap([2:n n:end-1],:));
%}


%% ERF_source_GA
clear all
erf_osc_datainfo;
datadir = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
k = 1;
    for subj = allsubs
        tmp{k} = load(fullfile([datadir, sprintf('sub-%03d_glm_parcel.mat', subj)]), 'tlck');
        k=k+1;
    end
    
    for k=1:32
        tlck{k} = tmp{k}.tlck;
    end
    clear tmp
    
    cfg=[];
    cfg.appenddim = 'rpt';
    tlck_all  = ft_appendtimelock(cfg, tlck{:});
   
    
    t1 = nearest(tlck_all.time, 0);
    t2 = nearest(tlck_all.time, 0.15);
    % align every parcel over subjects
   tlck_orig = tlck_all;
   for k = 1:length(tlck_all.label);
        erf = squeeze(tlck_all.trial(:,k,t1:t2));
        erf = diag(1./std(erf,[],2))*erf; % normalize before SVD
        
        [u,s,v]=svd(erf);
        u2(:,k) = sign(u(:,1));
    end
    u2 = repmat(u2, [1,1, length(tlck_all.time)]);
    
    tlck_all.trial = u2.*tlck_all.trial;
    
    
    cfg=[];
    cfg.channel = 'L_17_B05_03';
    L17 = ft_selectdata(cfg, tlck_all);    
    
    % ERF_all_source
    figure; plot(L17.time, squeeze(L17.trial));
    hold on;
    plot(L17.time, mean(squeeze(L17.trial)), 'LineWidth', 4);
    xlim([0 0.3])
    ylim([-3e-12 3e-12])
    
    % ERF_26_source
    figure; plot(L17.time, squeeze(L17.trial(25,1,:))); xlim([0 0.3])
%     cfg=[];
%     cfg.channel = {'L_17_B05_01', 'L_17_B05_02','L_17_B05_03', 'L_17_B05_04'};
%     v1_orig = ft_selectdata(cfg, tlck_orig);
%     v1_svd  = ft_selectdata(cfg, tlck_all);
%     cfg.channel = {'L_4_B05_03', 'L_4_B05_05','L_4_B05_08'}; % these look like they might correspond to hand
%     % otherwise use all 12 parcels
%     m1_orig = ft_selectdata(cfg, tlck_orig);
% %     m1_svd = ft_selectdata(cfg, tlck_all);
%     
%     cfgp.layout = 'layout_flatmap374.mat';
%     cfgp.gridscale=150;
%     cfgp.contournum=0;
% 
%     % v1
%     %ERF_source_m1_orig_GA 
% figure;
% hold on
% for i=1:32
% plot(v1_orig.time, squeeze(mean(v1_orig.trial(i,:,:),2)))
% end
% plot(v1_orig.time, squeeze(mean(mean(v1_orig.trial,2),1)), 'LineWidth', 4)
% %ERF_source_m1_svd_GA 
% figure;
% hold on
% for i=1:32
% plot(v1_svd.time, squeeze(mean(v1_svd.trial(i,:,:),2)))
% end
% plot(v1_svd.time, squeeze(mean(mean(v1_svd.trial,2),1)), 'LineWidth', 4)
% 
% % m1
% %ERF_source_m1_orig_GA 
% figure;
% hold on
% for i=1:32
% plot(m1_orig.time, squeeze(mean(m1_orig.trial(i,:,:),2)))
% end
% plot(m1_orig.time, squeeze(mean(mean(m1_orig.trial,2),1)), 'LineWidth', 4)
% %ERF_source_m1_svd_GA 
% figure;
% hold on
% for i=1:32
% plot(m1_svd.time, squeeze(mean(m1_svd.trial(i,:,:),2)))
% end
% plot(m1_svd.time, squeeze(mean(mean(m1_svd.trial,2),1)), 'LineWidth', 4)
% for subj=allsubs
%     tmp{subj} = load(sprintf('/project/3011085.02/results/erf/sub-%03d/erf_virtualchan_reversal.mat', subj);
%     tlck{subj} = ft_timelockanalysis([], tmp{subj}.tlck);
% end
% cfg=[];
% cfg.appenddim = 'rpt';
% tlck_all = ft_appendtimelock(cfg, tlck{allsubs});
% cfg=[];
% cfg.latency = [0 0.5];
% tlck_all = ft_selectdata(cfg, tlck_all);
% 
% plot(tlck_all.time, squeeze(tlck_all.trial))
% hold on
% plot(tlck_all.time, mean(squeeze(tlck_all.trial),1), 'LineWidth', 4)
% ylim([-0.02 0.02])
% 
% % ERF_source_2
% plot(tmp{2}.tlck.time, squeeze(tmp{2}.tlck.trial));
% hold on
% plot(tmp{2}.tlck.time, mean(squeeze(tmp{2}.tlck.trial),1), 'LineWidth', 4);
% ylim([-0.025 0.025])
% 
% % ERF_source_2
% plot(tmp{4}.tlck.time, squeeze(tmp{4}.tlck.trial));
% hold on
% plot(tmp{4}.tlck.time, mean(squeeze(tmp{4}.tlck.trial),1), 'LineWidth', 4);
% ylim([-0.03 0.03])


%% TFR_high_GA
clear all
erf_osc_datainfo;

% or load /project/3011085.02/results/TFR_all.mat
for subj=allsubs
tmp{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_high_onset.mat', subj));
end
for subj=allsubs
baseline{subj} = tmp{subj}.baseline;
tfa{subj} = tmp{subj}.tfa;
end
cfg=[];
cfg.appenddim='rpt';
tfa_all = ft_appendfreq(cfg, tfa{allsubs});
baseline_all = ft_appendfreq(cfg, baseline{allsubs});
baseline_all.powspctrm = repmat(mean(baseline_all.powspctrm,4),[1,1,1,71]);
baseline_all.time = tfa_all.time;
d=tfa_all;
d.powspctrm = tfa_all.powspctrm./baseline_all.powspctrm-1;



% % TFR_high_2
% cfgp=[];
% cfgp.layout='CTF275_helmet.mat';
% cfgp.xlim = [0 1.5];
% cfgp.zlim='maxabs';
% cfgp.colormap=flipud(brewermap(64,'RdBu'));
% cfgp.channel = {'MLO21', 'MLO22', 'MLO23', 'MLO31', 'MLO32', 'MLO33', 'MLO41', 'MLO42', 'MLO43', 'MLO52', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32', 'MRO41', 'MRO42', 'MRO43', 'MRO52', 'MZO02'};
% cfgp.trials = 2;
% ft_singleplotTFR(cfgp, d);
% 
% % TFR_topo_2
% cfgp=[];
% cfgp.layout='CTF275_helmet.mat';
% cfgp.xlim = [0.14 1.6];
% cfgp.ylim=[40 70];
% cfgp.zlim='maxabs';
% cfgp.marker = 'off';
% cfgp.highlight = 'on';
% cfgp.highlightchannel = {'MLO21', 'MLO22', 'MLO23', 'MLO31', 'MLO32', 'MLO33', 'MLO41', 'MLO42', 'MLO43', 'MLO52', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32', 'MRO41', 'MRO42', 'MRO43', 'MRO52', 'MZO02'};
% cfgp.highlightsymbol = 'O';
% cfgp.highlightsize = 8;
% cfgp.colormap=flipud(brewermap(64,'RdBu'));
% cfgp.trials = 2;
% ft_topoplotTFR(cfgp, d);

% TFR_26
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
% TFR_topo_26
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


% TFR_high_GA
cfgp=[];
cfgp.layout='CTF275_helmet.mat';
cfgp.xlim = [0 1.5];
cfgp.zlim=[-1 1];
n=11;
cmap = flipud(brewermap(2*n-1,'RdBu'));
cfgp.colormap=(cmap([2:n n:end-1],:));
cfgp.channel = {'MLO21', 'MLO22', 'MLO31', 'MRO21', 'MRO22', 'MRO31', 'MZO02'};
ft_singleplotTFR(cfgp, d);

% TFR_high_topo_GA
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

% TFA_high_tval_GA
% Nsub = length(allsubs);
% cfg             = [];
% cfg.method      = 'template'; % try 'distance' as well
% cfg.feedback    = 'no';
% neighbours      = ft_prepare_neighbours(cfg, tfa_all); % define neighbouring channels
% cfg                  = [];
% cfg.channel          = 'MEG';
% cfg.neighbours       = neighbours;
% cfg.parameter        = 'powspctrm';
% cfg.method           = 'analytic';
% cfg.statistic        = 'ft_statfun_depsamplesT';
% cfg.alpha            = 0.05;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2;
% stat = ft_freqstatistics(cfg, tfa_all, baseline_all);
% 
% 
% cfgp=[];
% cfgp.parameter = 'stat';
% cfgp.layout='CTF275_helmet.mat';
% cfgp.xlim = [0 1.5];
% cfgp.zlim='zeromax';
% cfgp.colormap=(brewermap(64,'OrRd'));
% cfgp.channel = {'MLO31', 'MLO22', 'MLO31', 'MRO22', 'MRO31', 'MZO02'};
% ft_singleplotTFR(cfgp, stat);

% TFA_high_tval_topo_GA
% cfgp=[];
% cfgp.parameter = 'stat';
% cfgp.layout='CTF275_helmet.mat';
% cfgp.xlim = [0.1 1.35];
% cfgp.ylim = [48 60];
% cfgp.zlim='maxabs';
% cfgp.colormap=flipud(brewermap(64,'RdBu'));
% ft_topoplotTFR(cfgp, stat);

%% TFR_high_source

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



%% Gamma boxplots
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

%% GLM channel
%{
clear all
erf_osc_datainfo;
load('/project/3011085.02/results/tlck_reversal_all.mat');
load('/project/3011085.02/results/GLM/jitter_covar/stat_glm_gamma_time_reversal.mat');

% stat_GLM_gamma_channel_topo_motor
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.zlim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.colorbar = 'yes';
cfgp.parameter = 'stat';
cfgp.xlim = [0.2333 0.3108];
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = find(any(stat.mask,2));
cfgp.highlightsymbol = '*';
figure; ft_topoplotER(cfgp, stat);
% and the correpsponding ERF:
% ERF_topo_motor_GA
cfgp = rmfield(cfgp, 'parameter');
figure; ft_topoplotER(cfgp, ft_timelockanalysis([],tlck_all));

% stat_GLM_gamma_channel_motor
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.ylim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.parameter = 'stat';
cfgp.channel = stat.label(find(any(stat.mask,2)));%{'MLC62', 'MLC63', 'MRC52', 'MRC54', 'MRC61', 'MRC62','MRC63','MZC01', 'MZC03', 'MZCC04'};
figure; ft_singleplotER(cfgp, stat);
vline(stat.time(find(any(stat.mask==1,1),1)))
vline(stat.time(find(any(stat.mask==1,1),1,'last')));
% and the correpsponding ERF:
% ERF_motor_GA
cfgp = rmfield(cfgp, 'parameter');
cfgp.xlim = [0 0.5];
figure; ft_singleplotER(cfgp, ft_timelockanalysis([],tlck_all));

% stat_GLM_gamma_channel_topo_p50
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.zlim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.colorbar = 'yes';
cfgp.parameter = 'stat';
cfgp.xlim = [0.035 0.05];
cfgp.marker = 'off';
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MLO32', 'MLO42', 'MLO43', 'MLO51', 'MLO52'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
figure;ft_topoplotER(cfgp, stat);
% and the correpsponding ERF:
% ERF_topo_p50_GA
cfgp.highlightchannel = {'MLO23', 'MLO32', 'MLO33', 'MLO43', 'MLO44'};
cfgp = rmfield(cfgp, 'parameter');
figure;ft_topoplotER(cfgp, ft_timelockanalysis([],tlck_all));

% stat_GLM_gamma_channel_p50
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.ylim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.parameter = 'stat';
cfgp.channel = {'MLO32', 'MLO42', 'MLO43', 'MLO51', 'MLO52'};
figure; ft_singleplotER(cfgp, stat);
vline(0.035)
vline(0.05)
% and the correpsponding ERF:
% ERF_p50_GA
cfgp = rmfield(cfgp, 'parameter');
cfgp.channel = {'MLO23', 'MLO32', 'MLO33', 'MLO43', 'MLO44'};
cfgp.xlim = [0 0.5];
figure; ft_singleplotER(cfgp, ft_timelockanalysis([],tlck_all));

% stat_GLM_gamma_channel_topo_p1
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.zlim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.colorbar = 'yes';
cfgp.parameter = 'stat';
cfgp.marker = 'off';
cfgp.xlim = [0.09 0.10];
cfgp.highlight = 'on';
cfgp.highlightchannel = {'MRC17', 'MRF46', 'MRF55', 'MRF56', 'MRF65', 'MRF67', 'MRP57', 'MRT11', 'MRT12', 'MRT13', 'MRT14'};
cfgp.highlightsymbol = 'O';
cfgp.highlightsize = 8;
ft_topoplotER(cfgp, stat);
% and the correpsponding ERF:
% ERF_topo_p1_GA
cfgp = rmfield(cfgp, 'parameter');
cfgp.highlightchannel = {'MLO32', 'MLO34', 'MLO42', 'MLO43', 'MLO44', 'MLT27', 'MLT37', 'MLT47', 'MLT57'};
figure; ft_topoplotER(cfgp, ft_timelockanalysis([],tlck_all));

% stat_GLM_gamma_channel_p1
cfgp=[];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.ylim = 'maxabs';
cfgp.colormap = flipud(brewermap(64, 'RdBu'));
cfgp.parameter = 'stat';
cfgp.channel = {'MRC17', 'MRF46', 'MRF55', 'MRF56', 'MRF65', 'MRF67', 'MRP57', 'MRT11', 'MRT12', 'MRT13', 'MRT14'};
ft_singleplotER(cfgp, stat);
vline(0.09);
vline(0.10);
% and the correpsponding ERF:
% ERF_p1_GA
cfgp.channel = {'MLO32', 'MLO34', 'MLO42', 'MLO43', 'MLO44', 'MLT27', 'MLT37', 'MLT47', 'MLT57'};
cfgp = rmfield(cfgp, 'parameter');
cfgp.xlim = [0 0.5];
figure; ft_singleplotER(cfgp, ft_timelockanalysis([],tlck_all));
%}
%% ERF source
clear
datadir1 = '/project/3011085.02/results/';
datadir2 = '/project/3011085.02/scripts/erfosc/analysis_JM_data/';
erf_osc_datainfo;
load(fullfile([datadir1, 'stat_glm_parcel.mat']), 'tlck_all', 'mu_erf');
k=1;
for subj=allsubs
    tmp{k} = load(fullfile([datadir2, sprintf('sub-%03d_parcel_blstd.mat', subj)]));
    tmp2{k} = tmp{k}.baseline_std;
    k=k+1;
end

baselinestd = permute(cat(2, tmp2{:}),[2,1]);

tmp = tlck_all;
tmp.trial = (tmp.trial-repmat(mu_erf, [1,1,360]))./repmat(baselinestd, [1,1,360]);

for k=1:32
    erf = squeeze(tlck_all.trial(k,:,:));
    [u,s,v] = svd(erf);
    u2 = permute(repmat(sign(u(:,1)), [1,1,360]),[2,1,3]);
    tmp.trial(k,:,:) = tmp.trial(k,:,:).*u2;
end

erf = squeeze(mean(tlck_all.trial,2));
[u,s,v] = svd(erf);
u2 = repmat(u(:,1), [1,370,360]);
tmp.trial = u2.*tmp.trial;

load('atlas_subparc374_8k.mat')
load cortex_inflated; atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); 

tlck_avg = [];
tlck_avg.label = atlas.parcellationlabel;
tlck_avg.brainordinate = atlas;
tlck_avg.pow = zeros(numel(tlck_avg.label),1);
tlck_avg.dimord = 'chan';

% ERF_source_1
tmp = tlck_all;
tmp.trial = (tmp.trial-repmat(mu_erf, [1,1,360]))./repmat(baselinestd, [1,1,360]);
t1 = nearest(tlck_all.time, 0.065);
t2 = nearest(tlck_all.time, 0.085);
for k=1:32
    erf = squeeze(tlck_all.trial(k,:,t1:t2));
    [u,s,v] = svd(erf);
    u2 = permute(repmat(sign(u(:,1)), [1,1,360]),[2,1,3]);
    tmp.trial(k,:,:) = tmp.trial(k,:,:).*u2;
end

erf = squeeze(mean(tlck_all.trial,2));
[u,s,v] = svd(erf);
u2 = repmat(u(:,1), [1,370,360]);
tmp.trial = u2.*tmp.trial;

load('atlas_subparc374_8k.mat')
load cortex_inflated; atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); 

tlck_avg = [];
tlck_avg.label = atlas.parcellationlabel;
tlck_avg.brainordinate = atlas;
tlck_avg.pow = zeros(numel(tlck_avg.label),1);
tlck_avg.dimord = 'chan';


tlck_avg.pow(selparc,1) = mean((mean(tmp.trial(:,:,t1:t2),3)),1);

cfg=[];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap = brewermap(64, 'OrRd');%flipud(brewermap(64,'RdBu'));
cfg.colorbar = 'yes';
cfg.maskstyle = 'colormix';
cfg.maskparameter = 'pow';
cfg.camlight = 'no';
cfg.funcolorlim = [0 0.03];%[-0.025 0.025];%'maxabs';
ft_sourceplot(cfg, tlck_avg);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([40 6])

% ERF_source_2
tmp = tlck_all;
tmp.trial = (tmp.trial-repmat(mu_erf, [1,1,360]))./repmat(baselinestd, [1,1,360]);
t1 = nearest(tlck_all.time, 0.12);
t2 = nearest(tlck_all.time, 0.16);
for k=1:32
    erf = squeeze(tlck_all.trial(k,:,t1:t2));
    [u,s,v] = svd(erf);
    u2 = permute(repmat(sign(u(:,1)), [1,1,360]),[2,1,3]);
    tmp.trial(k,:,:) = tmp.trial(k,:,:).*u2;
end

erf = squeeze(mean(tlck_all.trial,2));
[u,s,v] = svd(erf);
u2 = repmat(u(:,1), [1,370,360]);
tmp.trial = u2.*tmp.trial;

load('atlas_subparc374_8k.mat')
load cortex_inflated; atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); 

tlck_avg = [];
tlck_avg.label = atlas.parcellationlabel;
tlck_avg.brainordinate = atlas;
tlck_avg.pow = zeros(numel(tlck_avg.label),1);
tlck_avg.dimord = 'chan';

tlck_avg.pow(selparc,1) = mean((mean(tmp.trial(:,:,t1:t2),3)),1);

cfg=[];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap = brewermap(64, 'OrRd');%flipud(brewermap(64,'RdBu'));
cfg.colorbar = 'yes';
cfg.maskstyle = 'colormix';
cfg.maskparameter = 'pow';
cfg.camlight = 'no';
cfg.funcolorlim = [0 0.03];%[-0.025 0.025];%'maxabs';
ft_sourceplot(cfg, tlck_avg);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([62 16])

% ERF_source_3
%{
% t1 = nearest(tlck_all.time, 0.2);
% t2 = nearest(tlck_all.time, 0.24);
% tmp = tlck_all;
% tmp.trial = (tmp.trial-repmat(mu_erf, [1,1,360]))./repmat(baselinestd, [1,1,360]);
% 
% for k=1:32
%     erf = squeeze(tlck_all.trial(k,:,t1:t2));
%     [u,s,v] = svd(erf);
%     u2 = permute(repmat(sign(u(:,1)), [1,1,360]),[2,1,3]);
%     tmp.trial(k,:,:) = tmp.trial(k,:,:).*u2;
% end
% 
% erf = squeeze(mean(tlck_all.trial,2));
% [u,s,v] = svd(erf);
% u2 = repmat(u(:,1), [1,370,360]);
% tmp.trial = u2.*tmp.trial;
% 
% load('atlas_subparc374_8k.mat')
% load cortex_inflated; atlas.pos=ctx.pos;
% exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
% selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); 
% 
% tlck_avg = [];
% tlck_avg.label = atlas.parcellationlabel;
% tlck_avg.brainordinate = atlas;
% tlck_avg.pow = zeros(numel(tlck_avg.label),1);
% tlck_avg.dimord = 'chan';
% tlck_avg.pow(selparc,1) = mean((mean(tmp.trial(:,:,t1:t2),3)),1);
% 
% cfg=[];
% cfg.method = 'surface';
% cfg.funparameter = 'pow';
% cfg.funcolormap = flipud(brewermap(64,'RdBu'));
% cfg.colorbar = 'yes';
% cfg.maskstyle = 'colormix';
% cfg.maskparameter = 'pow';
% cfg.camlight = 'no';
% cfg.funcolorlim = [-0.025 0.025];%'maxabs';
% ft_sourceplot(cfg, tlck_avg);
% h = light('position', [-1 0 -0.1]);
% h2=light; set(h2, 'position', [1 0 -0.1]);
% material dull;view([68 12])

% ERF_source_4
tmp = tlck_all;
t1 = nearest(tlck_all.time, 0.25);
t2 = nearest(tlck_all.time, 0.35);
tmp.trial = (tmp.trial-repmat(mu_erf, [1,1,360]))./repmat(baselinestd, [1,1,360]);
for k=1:32
    erf = squeeze(tlck_all.trial(k,:,t1:t2));
    [u,s,v] = svd(erf);
    u2 = permute(repmat(sign(u(:,1)), [1,1,360]),[2,1,3]);
    tmp.trial(k,:,:) = tmp.trial(k,:,:).*u2;
end

erf = squeeze(mean(tlck_all.trial,2));
[u,s,v] = svd(erf);
u2 = repmat(u(:,1), [1,370,360]);
tmp.trial = u2.*tmp.trial;

load('atlas_subparc374_8k.mat')
load cortex_inflated; atlas.pos=ctx.pos;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); 

tlck_avg = [];
tlck_avg.label = atlas.parcellationlabel;
tlck_avg.brainordinate = atlas;
tlck_avg.pow = zeros(numel(tlck_avg.label),1);
tlck_avg.dimord = 'chan';
tlck_avg.pow(selparc,1) = -mean((mean(tmp.trial(:,:,t1:t2),3)),1);


cfg=[];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap = brewermap(64, 'OrRd');%flipud(brewermap(64,'RdBu'));
cfg.colorbar = 'yes';
cfg.maskstyle = 'colormix';
cfg.maskparameter = 'pow';
cfg.camlight = 'no';
cfg.funcolorlim = [0 0.03];%[-0.025 0.025];%'maxabs';
ft_sourceplot(cfg, tlck_avg);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([-71 20])
%}

%% Correlation Gamma ERF
load('/project/3011085.02/results/stat_peakpicking3.mat', 'stat','S');
load atlas_subparc374_8k.mat
load cortex_inflated_shifted; atlas.pos=ctx.pos;
stat.brainordinate=atlas;


% required to hack around a little bit. For some reason the anatomical data
% is not plotted as should. Workaround: plot atlas first, and onto this the
% functional data (comment out line line 705 in ft_sourceplot, which opens
% a new figure).

cfgx = [];
cfgx.method='surface';
cfgx.funparameter='stat';
cfgx.funcolormap = flipud(brewermap(64, 'RdBu'));

% n=7;
% cmap = flipud(brewermap(2*n-1,'RdBu'));
% cfgx.funcolormap=(cmap([2:n n:end-1],:));

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
x=find(stat.posclusterslabelmat==1);
y = mean(S.rho(:,x),2);
scatter(1:32, sort(y));
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
y = mean(rho(:,x),2);
scatter(1:32, sort(y));
addpath /project/3011085.02/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
figure; iosr.statistics.boxPlot(y)

%% GLM source parcels
%{
clear
savedir = '/project/3011085.02/results/';
load(fullfile([savedir, 'stat_glm_parcel.mat']), 'stat','betas_align_norm', 'tlck_all');

load('atlas_subparc374_8k.mat')
load cortex_inflated_shifted; atlas.pos=ctx.pos;

beta_stat = [];
beta_stat.label = atlas.parcellationlabel;
beta_stat.brainordinate = atlas;
beta_stat.pow = zeros(numel(beta_stat.label),1);
beta_stat.dimord = 'chan';

exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label);

beta_stat.pow(selparc,1) = mean(stat.stat.*stat.mask,2);

% stat_glm_parc_1a
cfg=[];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap = (brewermap(64,'OrRd'));
cfg.colorbar = 'no';
cfg.maskstyle = 'colormix';
cfg.maskparameter = 'pow';
cfg.camlight = 'no';
cfg.funcolorlim = 'zeromax';
ft_sourceplot(cfg, beta_stat);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([91 18])

% stat_glm_parc_1b
ft_sourceplot(cfg, beta_stat);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([-96 -2])

% stat_glm_parc_1c
figure;
plot(stat.time, mean(stat.mask.*stat.stat), 'LineWidth', 3)


%% GLM source parcels - response locked
clear
savedir = '/project/3011085.02/results/';
load(fullfile([savedir, 'stat_glm_parcelresp.mat']));

load('atlas_subparc374_8k.mat')
load cortex_inflated_shifted; atlas.pos=ctx.pos;

beta_stat = [];
beta_stat.label = atlas.parcellationlabel;
beta_stat.brainordinate = atlas;
beta_stat.pow = zeros(numel(beta_stat.label),1);
beta_stat.dimord = 'chan';

exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label);

beta_stat.pow(selparc,1) = nanmean(stat.stat.*stat.mask,2);

% stat_glm_parc_1a
cfg=[];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.funcolormap = (brewermap(64,'OrRd'));
cfg.colorbar = 'yes';
cfg.maskstyle = 'colormix';
cfg.maskparameter = 'pow';
cfg.camlight = 'no';
cfg.funcolorlim = 'zeromax';
ft_sourceplot(cfg, beta_stat);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([91 18])

% stat_glm_parc_1b
ft_sourceplot(cfg, beta_stat);
h = light('position', [-1 0 -0.1]);
h2=light; set(h2, 'position', [1 0 -0.1]);
material dull;view([-96 -2])

% stat_glm_parc_1c
figure;
plot(stat.time, mean(stat.mask.*stat.stat), 'LineWidth', 3)




% sourcemovie
% load cortex_inflated.mat
% load('atlas_subparc374_8k.mat')
% load cortex_inflated_shifted;
% source = tlck_svd2;
% source.brainordinate = atlas;
% source.brainordinate.parcellationlabel = source.label;
% source.brainordinate.pos = ctx.pos;
% 
% cfg=[];
% cfg.funparameter = 'trial';
% cfg.funcolormap = flipud(brewermap(64, 'RdBu'));
% cfg.parcellation = 'parcellation';
% ft_sourcemovie(cfg, source);
%}
