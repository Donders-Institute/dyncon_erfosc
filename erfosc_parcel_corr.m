clear
erf_osc_datainfo;
datadir = '/project/3011085.02/analysis/';

k=1;
for subj=allsubs
    tmp{k} = load(fullfile([datadir, 'erf/', sprintf('sub-%03d/sub-%03d_parcel_tlck', subj, subj)]), 'erfdata');
    tlck{k} = tmp{k}.erfdata;
    tmp2{k} = load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel.mat', subj, subj), 'gammaPow');
    gammaPow{k} = tmp2{k}.gammaPow;
    k=k+1;
end

for k=1:32
    cfg=[];
    cfg.channel = {'L_17*', 'R_17*'};
    cfg.latency = [0 0.15];
    tlck{k} = ft_selectdata(cfg, tlck{k});
end
t1 = nearest(tlck{1}.time, 0.06);

for k=1:32
    erf = squeeze(mean(tlck{k}.trial(:,:,t1:end),1));
    [u,s,v]=svd(erf);
    u2 = sign(u(:,1));
    u2 = permute(repmat(u2, [1 size(tlck{k}.trial,1) length(tlck{k}.time)]), [2,1,3]);
    tlck{k}.trial = u2.*tlck{k}.trial;
    clear u2 u s v erf
end

for k=1:32
    erf{k} = squeeze(mean(tlck{k}.trial,2));
end

peak = [0.09208, 0.08875, 0.1121, 0.1004, 0.07042, 0.1304, 0.1021, 0.07708, 0.09708, 0.07042, 0.09875, 0.08875, 0.08875, 0.08208, 0.07542, 0.1254, 0.08375, ...
    0.09708, 0.07375, 0.08708, 0.07708, 0.07208, 0.07042, 0.08042, 0.07542, 0.09208, 0.1221, 0.06542, 0.06542, 0.06542, 0.08375, 0.06542];

for k=1:32
idx(k) = nearest(tlck{1}.time, peak(k));
end

for k=1:32
    amp{k} = mean(erf{k}(:,idx(k)-3:idx(k)+3),2);
end

for k=1:32
[r(k), p(k)] = corr(amp{k}, gammaPow{k}', 'type', 'Spearman');
end

effect = ttest(r)

scatter(1:32, r);









