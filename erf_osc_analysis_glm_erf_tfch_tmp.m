function erf_osc_analysis_glm_erf_tfch(subj, isPilot, freqRange, zeropoint, erfoi, doDSS)
% linear regression of peak amplitude over time-frequency (with fixed
% channel) or over frequency-channel (with fixed (avg) time).

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end
if nargin<3 || isempty(freqRange)
    freqRange = 'high'; % can be 'high' or 'low'; depending on which frequency range you want to do regression (2-30 Hz or 28-100 Hz)
end
if nargin<4 || isempty(zeropoint)
    zeropoint = 'onset'; % other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if nargin<5 || isempty(erfoi) % erf of interest
    erfoi = 'reversal'; % other option 'onset': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if nargin<6 || isempty(doDSS)
  doDSS=false;
end

% Initiate Diary
%ft_diary('on')

%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    %load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    %load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s.mat', subj, zeropoint));
    if strcmp(erfoi, 'reversal')
        if doDSS
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'reversal', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    else
        if doDSS
            [data, nComp_keep] = erf_osc_analysis_dss(subj,isPilot, 'onset', false);
        else
            load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/cleandata.mat', subj));
        end
    end
end
%% if no data_dss cleaning is done beforehand, do this
if ~doDSS
    data=dataClean;
    idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
    nTrials = length(idxM);
    
    cfg=[];
    cfg.trials = idxM;
    cfg.channel = 'MEG';
    data = ft_selectdata(cfg, data);
    
    if strcmp(erfoi, 'reversal')
        cfg=[];
        cfg.offset = -(data.trialinfo(:,5)-data.trialinfo(:,4));
        data = ft_redefinetrial(cfg, data);
    end
    clear dataClean;
else
    data=data_dss;
    clear data_dss;
end

fs=data.fsample;
nTrials = length(data.trial);


%% select p1 window for regression
cfg=[];
cfg.vartrllength=2;
tlck = ft_timelockanalysis(cfg, data);

t1p1 = nearest(tlck.time, 0.06);
t2p1 = nearest(tlck.time, 0.12);
cfg=[];
cfg.channel = {'MRO', 'MRP', 'MLO', 'MRO', 'MZO', 'MZP'};
cfg.latency = [tlck.time(t1p1) tlck.time(t2p1)];
tlck = ft_selectdata(cfg, tlck);

time = tlck.time;
[~, maxchan] = max(abs(mean(tlck.avg,2))); % find channel with max amplitude
% calculate mean over every window to find out which window has the maximum
% mean amplitude (for the maximum channel!). Take the mean amplitude in
% this latency window as regression weight.
maxchanid = tlck.label(maxchan);
maxchanidx = find(strcmp(maxchanid, data.label));
halfwindowlength = 8;
i=1;
for ntrl = halfwindowlength+1:length(time)-halfwindowlength
    win(i,:) = [time(ntrl-8), time(ntrl+8)];
    avg(i,1) = mean(tlck.avg(maxchan,ntrl-8:ntrl+8),2);
    i=i+1;
end
[~, window] = max(abs(avg));
lat = win(window,:);


%% Regression p1 amplitude over time-frequency-channel
cfg=[];
cfg.latency = lat;%[-1.5 0.65]; -> JM: dit bespaart je een hoop geheugen, toch? Zeker als je zo een cat(3, ...) doet
cfg.avgovertime = 'yes';
data = ft_selectdata(cfg, data);
p1amp = cat(2,data.trial{:});

%cfg=[];
%cfg.latency = [-1.5 0.65];
%trialdata = cat(3,data.trial{:});

%idxtime = [nearest(data.time{1}, lat(1)) : nearest(data.time{1}, lat(2))];
%p1amp = squeeze(mean(trialdata(:,idxtime,:),2));

% baselinecorrect with average baseline over trials
if isPilot
  load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
  load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%s.mat', subj, zeropoint));
end
if strcmp(freqRange, 'high')
  tfa = tfaHigh;
  clear tfaHigh tfaLow; % JM this step only renames the variables, but I think this should be done.
  % it prevents code from being duplicated in the subsequent step, just do
  % to a difference in variable name. for the next time I would solve this
  % by storing the low and high ranges in a separate file. Note: in the
  % current situation you could consider to read selectively the variable
  % of interest, (as per the loading line above load(filename, 'tfaHigh');
  % or so. 
elseif strcmp(freqRange, 'low')
  tfa = tfaLow;
  clear tfaLow tfaHigh;
else
  %throw an informative error
  keyboard
end


if strcmp(zeropoint, 'onset')
    cfg=[];
    cfg.latency = [-1 1]; % shortest baseline window is 1 second, shortest reversal 1 second as well
    tfa = ft_selectdata(cfg, tfa);
else
    cfg=[];
    cfg.latency = [-1.5 0];
    tfa = ft_selectdata(cfg, tfa);
end

nchan = size(tfa.powspctrm, 2);
nfreq = size(tfa.powspctrm, 3);
ntime = size(tfa.powspctrm, 4);
ntrl = size(p1amp,2);
for freq=1:nfreq
    for ch=1:nchan
        %JM: hier ben ik toch niet zo blij mee. D.w.z. de kanaalspecifieke
        %p1amp. Voor kanalen zonder een duidelijke ERF zal dit gewoon ruis
        %worden, en daardoor niet informatief (en logisch dat de betas rond
        %0 gaan zijn. Ik wil graag dat je een versie van de resultaten
        %maakt, waarbij je als p1amp een verdedigbaar groepje kanalen
        %averaged, waarbij het groepje kanalen genomen is uit de subset die
        %je boven gebruikt voor de latency bepaling. (bijvoorbeeld het
        %gemiddelde van de top 5%).
        design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
        design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
        Y_h = squeeze(squeeze(tfa.powspctrm(:,ch,freq,:)));
        betas(freq,ch,:,:) = design'\Y_h;
    end
end

design_bias = ones(1, ntrl);
p1amp_norm = (p1amp - repmat(mean(p1amp,2), [1 ntrl]))./repmat(std(p1amp, [], 2), [1, ntrl]);

numShuffles=100;
all_shuffles = zeros(nchan, nfreq, ntime, numShuffles);

% prepare a data object and cfg for the planar transform
tl=[];
tl.avg    = squeeze(betas(:,:,2,:));
tl.time   = tfa.time;
tl.dimord = 'rpt_chan_time';
tl.label  = tfa.label;
tl.grad   = tfa.grad;

cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
cfg.planarmethod    = 'sincos';


for iShuffle = 1:numShuffles
  shufvec = randperm(ntrl); % JM:probably it's best to do the randomization the same
  % for all channels/times/frequencies, this requires the shuffles to be in
  % the outer loop
  
  betas_shuf = zeros(nchan, nfreq*ntime); % store betas for this shuffle across channels
  for ch = 1:nchan
    design_shuf = [design_bias; p1amp_norm(ch,:)];
    design_shuf = design_shuf(:,shufvec);
    Y1  = reshape(tfa.powspctrm(:,ch,:,:),[],nfreq*ntime);
    tmp = design_shuf'\Y1;
    betas_shuf(ch,:) = tmp(2,:);
  end
  tl.time = 1:nfreq*ntime;
  tl.avg  = betas_shuf;
  tmp     = ft_combineplanar([], ft_megplanar(cfg, tl));
  betas_shuf = tmp.avg;
  all_shuffles(:,:,:,iShuffle) = reshape(betas_shuf, [nchan nfreq ntime]);
end
avg_shuffles = mean(all_shuffles,4);
std_shuffles = std(all_shuffles,[],4);

% JM weet niet zeker of the all_shuffles ook ge saved moeten worden, nu de
% avg er mogelijk wat gezonder uit gaat zien.

%% planar gradient transformation of beta weights
% make timelocked structure where planar gradient transformation can be
% applied to (that's why dimord is strange)
tl=[];
tl.avg    = squeeze(betas(:,:,2,:));
tl.time   = tfa.time;
tl.dimord = 'rpt_chan_time';
tl.label  = tfa.label;
tl.grad   = tfa.grad;

% also put betas in a time-freq structure
tf = rmfield(tl,'avg');
tf.dimord = 'chan_freq_time';
tf.powspctrm = permute(tl.avg, [2,1,3]);
tf.freq = tfa.freq;

% planar combination
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
cfg.planarmethod    = 'sincos';
tlPlanar            = ft_megplanar(cfg, tl);

cfg           = [];
bPlanarCmb  = ft_combineplanar(cfg,tlPlanar);

% put it back in a freq-data structure
bPlanarCmb.powspctrm = permute(bPlanarCmb.trial, [2,1,3]);
bPlanarCmb           = rmfield(bPlanarCmb, 'trial');
bPlanarCmb.freq      = tfa.freq;
bPlanarCmb.dimord    = 'chan_freq_time';

tl_shuffle_avg      = rmfield(tl, 'avg');
tl_shuffle_avg.avg  = permute(avg_shuffles, [2,1,3]);
tl_shuffle_std      = rmfield(tl, 'avg');
tl_shuffle_std.avg  = permute(std_shuffles,[2,1,3]);

cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl_shuffle_avg);
cfg.planarmethod    = 'sincos';
tl_shuffle_avgPlanar= ft_megplanar(cfg, tl_shuffle_avg);
tl_shuffle_stdPlanar= ft_megplanar(cfg, tl_shuffle_std);

cfg                      = [];
tl_shuffle_avgCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_avgPlanar);
tl_shuffle_stdCmbPlanar  = ft_combineplanar(cfg,tl_shuffle_stdPlanar);

shuffle_avgCmbPl = rmfield(bPlanarCmb, 'powspctrm');
shuffle_avgCmbPl.powspctrm = permute(tl_shuffle_avgCmbPlanar.trial, [2,1,3]);
shuffle_stdCmbPl = rmfield(bPlanarCmb, 'powspctrm');
shuffle_stdCmbPl.powspctrm = permute(tl_shuffle_stdCmbPlanar.trial, [2,1,3]);

%% Save

if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_tf_%s_%s_erf_%s', subj, freqRange, zeropoint, erfoi);
else
    filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tf_%s_%s_erf_%s', subj, freqRange, zeropoint, erfoi);
end

save(fullfile([filename '.mat']), 'betas','bPlanarCmb','tf','lat','maxchanid','maxchanidx', 'shuffle_avgCmbPl', 'shuffle_stdCmbPl', 'p1amp', '-v7.3');



ft_diary('off')

