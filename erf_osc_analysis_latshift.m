% for each parcel, divide the trials in bins according to reaction time.
% Then, for each parcel, convolve the time series of the bin average with
% that of the next bin. The time at which this result peaks, is the delay
% between the maxima of the two bins (correct?? or should this be divided
% by two?). find the maxima of the convolved signal and correlate this with
% the mean reaction time per bin.
for subj=[1:9 11:15]
load(sprintf('/project/3011085.02/scripts/erfosc/analysis_JM_data/sub-%03d_erfparc.mat', subj))
load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj));

fs = parceldata_shift.fsample;

% order the trials according to trial number
[tmp, idx] = sort(parceldata_shift.trialinfo);
parceldata_shift.trialinfo = parceldata_shift.trialinfo(idx);
parceldata_shift.trial = parceldata_shift.trial(idx);
clear idx

% timelock analysis
cfg=[];
cfg.keeptrials= 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
cfg.preproc.demean = 'yes';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 20;
cfg.preproc.lpfilttype = 'firws';
tlck = ft_timelockanalysis(cfg, parceldata_shift);

% ntrials = numel(parceldata_shift.trial);
% % find the parcels in which there is activation
% cfg=[];
% cfg.latency = [0 0.5];
% t1  =ft_selectdata(cfg, tlck);
% cfg=[];
% cfg.latency = [-0.2 0];
% cfg.avgovertime = 'yes';
% t2  =ft_selectdata(cfg, tlck);
% t2.trial = repmat(t2.trial, [1 1 300]);
% t2.time = t1.time;
% cfg=[];
% cfg.method = 'montecarlo';
% cfg.statistic = 'depsamplesT';
% cfg.uvar = 2;
% cfg.ivar=1;
% cfg.design = [ones(1,ntrials), 2*ones(1,ntrials); 1:ntrials, 1:ntrials];
% cfg.alpha = 0.05;
% cfg.correctm = 'no';
% cfg.numrandomization  =1000;
% stat = ft_timelockstatistics(cfg, t1, t2);




% devide the RT distribution in 9 bins according to reaction time
u = mean(rt);
sd = std(rt);
nbins=7;
[Y, edges] = discretize([u-2*sd u+2*sd], nbins); % take all RT within 2 SD of the mean

X = find(rt>=u-2*sd & rt<=u+2*sd);
rtsel = rt(X);
nbins_orig=nbins;
whichbin=[];
for k=1:nbins
    if nbins~=nbins_orig && k==nbins_orig
        break
    end
    idx = find(rtsel>=edges(k) & rtsel<edges(k+1));
    if isempty(idx)
        edges(k)=[];
        nbins=nbins-1;
    end
    clear idx
end

cfg=[];
cfg.trials = X;
tlck = ft_selectdata(cfg, tlck);

% select trials belonging to each bin
cfg=[];
cfg.latency = [0 0.7];
cfg.avgoverrpt = 'yes';
for k=1:nbins
    idx = find(rtsel>=edges(k) & rtsel<edges(k+1));
    rtbin{k} = rtsel(idx);
    cfg.trials = idx;
    tlckbin{k} = ft_selectdata(cfg, tlck);
end

% mean RT for each bin
for k=1:nbins
    rtbinavg(k,1) = mean(rtbin{k});
end

% difference between each bin's avg and the avg in bin 1 (difference in
% response time between fastest trials and slower trials)
for k=1:nbins-1
    rtdelay(k) = rtbinavg(k+1)-rtbinavg(1);
end

time1 = tlckbin{1}.time;
time2 = sort(unique([-tlckbin{1}.time(2:end) tlckbin{1}.time]));
for j=1:numel(tlck.label)
    s_j = zeros(nbins, numel(time1));
    % select parcel
    for k=1:nbins
        s_j(k,:) = tlckbin{k}.trial(j,:);
    end
    s_j_perm = s_j(randperm(size(s_j,1)),:);
    % convolve the time course of a binavg with the next binavg
    for k=1:nbins-1
        binconv(k,:) = conv(s_j(k,:),  s_j(k+1,:));
        binconv_perm(k,:) = conv(s_j_perm(k,:),  s_j_perm(k+1,:));
    end
    binconv = binconv(:, (size(binconv,2)+1)/2:end); % only take the
    % positive values because we're interested in the DELAY of slower
    % trials wrt the fastest trials
    binconv_perm = binconv_perm(:, (size(binconv_perm,2)+1)/2:end);
    
    [val maxidx] = max(binconv');
    ampdelay{j} = time1(maxidx);
    
    [val2 maxidx2] = max(binconv_perm');
    ampdelay_perm{j} = time1(maxidx2);
    
    [rho(j) pval(j)] = corr(rtdelay', ampdelay{j}');
    [rho_perm(j) pval_perm(j)] = corr(rtdelay', ampdelay_perm{j}');
    clear s_j dum binconv binconv_perm s_j_perm
end

% convolve the middle bin for each parcel with that of the mean

cfg=[];
cfg.latency = [0 0.7];
avgbin = ft_selectdata(cfg,tlck);
avgbin.trial = squeeze(mean(avgbin.trial));%tlckbin{round(nbins/2)}; % take the middle bin time course
avgbin.trial = abs(avgbin.trial); % get rid of polarity ambiguity
parcelconv = zeros(numel(avgbin.label), numel(time2));
idxl = match_str(avgbin.label, 'L_17_B05_01');
idxl = match_str(avgbin.label, 'R_17_B05_01');

for k=1:numel(avgbin.label)
    if match_str(avgbin.label(k), 'L_17_B05_01') | match_str(avgbin.label(k), 'R_17_B05_01')
        % do nothing
    elseif match_str(avgbin.label{1}(1), 'L') % compare with left hemisphere V1 parcel
        parcelconv(k,:) = conv(avgbin.trial(idxl,:), avgbin.trial(k,:));
    elseif match_str(avgbin.label{1}(1), 'R') % compare with right hemisphere V1 parcel
        parcelconv(k,:) = conv(avgbin.trial(idxr,:), avgbin.trial(k,:));
    end
end
parcelconv = parcelconv(:, (size(parcelconv,2)+1)/2:end); % take only positive delays
[val maxidx] = max(parcelconv');
parceldelay = time1(maxidx);

S.rho = rho;
load atlas_subparc374_8k.mat
load cortex_inflated_shifted.mat;
atlas.pos=ctx.pos;
S.label = atlas.parcellationlabel;
S.dimord = 'chan_freq';
S.freq=1;
S.brainordinate = atlas;
S.rho =  zeros(374,1);
S_perm = S;
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
S.rho(setdiff(1:374, exclude_label)) = rho;
S.pval =  zeros(374,1);
S.pval(setdiff(1:374, exclude_label)) = pval;
S.parceldelay = zeros(374,1);
S.parceldelay(setdiff(1:374, exclude_label)) = parceldelay;

S_perm.rho = zeros(374,1);
S_perm.rho(setdiff(1:374, exclude_label)) = rho_perm;
S_perm.pval =  zeros(374,1);
S_perm.pval(setdiff(1:374, exclude_label)) = pval_perm;
% cfgp.funparameter = 'rho';
% cfgp.method = 'surface';
% cfgp.funcolorlim = 'zeromax';
% ft_sourceplot(cfgp, S);


filename = sprintf('/project/3011085.02/scripts/erfosc/analysis_JM_data/sub-%03d_latshift.mat', subj);
save(filename, 'S', 'rtbin', 'S_perm');
keep subj
end