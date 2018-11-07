
erf_osc_datainfo;
erfoi='reversal';

for subj=allsubs

g{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/allori/gamma_virtual_channel.mat', subj), 'gammaPow');
g{subj} = g{subj}.gammaPow;
data{subj} = erf_osc_analysis_lcmv_orientation(subj, erfoi); 
end


%% define latency
cfg=[];
cfg.latency = [-0.1 0.5];
cfg2=[];
cfg2.vartrllength = 2;
cfg2.preproc.baselinewindow = [-0.1 0];
cfg2.preproc.demean = 'yes';
cfg2.keeptrials = 'yes';

for subj=allsubs
        data{subj} = ft_selectdata(cfg, data{subj});
        tlck{subj} = ft_timelockanalysis(cfg2, data{subj});
end
time = tlck{1}.time;
a = nearest(time, 0);
b = nearest(time, 0.3);
% latency of highest peak in [0 300ms]
k=1;
for subj=allsubs
latency(k) = max(mean(abs(squeeze(tlck{subj}.trial(:,1,a:b))),1));
k=k+1;
end
% latency of first peak
% latency = [0.0925 0.08 0.07833 0.07083 0.07425 0.07333 0.0725 0.0725 0.0858 0.07 0.1008 0.0725 0.0625 0.1325 0.06583 0.06917 0.1392 0.04667 0.04417 0.1083 0.0875 0.045 0.425 0.04917 0.07 0.09167 0.04917 0.09333 0.1142 0.07083 0.07667 0.07333];

k=1;
% amplitude is average amplitude within [peak-5/fs peak+5/fs]
for subj=allsubs
    t1(k) = nearest(time, latency(k)) - 5;
    t2(k) = nearest(time, latency(k)) + 5;
    amp{subj}(:,1) = mean(tlck{subj}.trial(:,1,t1:t2),3);
    k = k+1;
end

k=1;
for subj=allsubs
    R(k) = corr(g{subj}', amp{subj}, 'type', 'Spearman');
    k = k+1;
end


%% dipole flip
for subj=allsubs
    tlck{subj} = rmfield(tlck{subj},'trial');
end

cfg               = [];
cfg.appenddim     = 'rpt';
tlck_GA = ft_appendtimelock(cfg, tlck{allsubs});

cfg=[];
cfg.latency = [0 0.5];
tlck_GA = ft_selectdata(cfg, tlck_GA);

% dipole flip
erf = squeeze(tlck_GA.trial);
erf = diag(1./std(erf,[],2))*erf; % normalize before SVD

[u,s,v]=svd(erf);

u2 = sign(u(:,1))';
R_flipped = u2.*R;




