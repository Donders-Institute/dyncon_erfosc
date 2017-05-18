function erf_osc_analysis_glm_erf_tfch(subj, isPilot, freqRange, zeropoint;)
% linear regression of peak amplitude over time-frequency (with fixed
% channel) or over frequency-channel (with fixed (avg) time).

if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = false;
end
if isempty(isPilot);
    isPilot = false;
end
if nargin<3
    freqRange = 'high'; % can be 'high' or 'low'; depending on which frequency range you want to do regression (2-30 Hz or 28-100 Hz)
end
if isempty(freqRange)
    freqRange = 'high';
end
if nargin<4
    zeropoint = 'onset'; % other option 'reversal': redefine time axis to stimulus reversal or keep it at stimulus onset
end
if isempty(zeropoint);
    zeropoint = 'onset';
end

%% Initiate Diary
workSpace = whos;
diaryname = tempname(fullfile([getenv('HOME'), '/tmp']));
diary(diaryname) % save command window output
fname = mfilename('fullpath')
datetime

fid = fopen(fullfile([fname '.m']));
tline = fgets(fid); % returns first line of fid
while ischar(tline) % at the end of the script tline=-1
    disp(tline) % display tline
    tline = fgets(fid); % returns the next line of fid
end
fclose(fid);

for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
end


%% load data
erf_osc_datainfo;
if isPilot
    load(sprintf('/project/3011085.02/results/erf/pilot-%03d/dss.mat', subj), 'data_dss');
    load(sprintf('/project/3011085.02/results/freq/pilot-%03d/gamma_virtual_channel.mat', subj), 'gammaChan');
else
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa_%d.mat', subj, zeropoint));
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/dss.mat', subj), 'data_dss');
    %     load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
end
fs=data_dss.fsample;
nTrials = length(data_dss.trial);

%% select p1 window for regression

tlck = ft_timelockanalysis([], data_dss);
t1p1 = nearest(tlck.time, 0.06);
t2p1 = nearest(tlck.time, 0.12);
cfg=[];
cfg.latency = [tlck.time(t1p1) tlck.time(t2p1)];
tlck = ft_selectdata(cfg, tlck);


time = tlck.time;
[~, maxchan] = max(abs(mean(tlck.avg,2))); % find channel with max amplitude
% calculate mean over every window to find out which window has the maximum
% mean amplitude (for the maximum channel!). Take the mean amplitude in
% this latency window as regression weight.
halfwindowlength = 8;
i=1;
for t = halfwindowlength+1:length(time)-halfwindowlength;
    win(i,:) = [time(t-8), time(t+8)];
    avg(i,1) = mean(tlck.avg(maxchan,t-8:t+8),2);
    i=i+1;
end
[~, window] = max(abs(avg));
lat = win(window,:);

%% Regression p1 amplitude over time-frequency-channel

trialdata = cat(3,data_dss.trial{:});

idxtime = [nearest(data_dss.time{1}, lat(1)) : nearest(data_dss.time{1}, lat(2))];
p1amp = squeeze(mean(trialdata(:,idxtime,:),2));

% baselinecorrect with average baseline over trials
if strcmp(freqRange, 'high');
    cfg=[];
    cfg.latency = [-1+1/fs -0.25];
    cfg.avgoverrpt = 'yes';
    cfg.avgovertime = 'yes';
    baselineH = ft_selectdata(cfg, tfaHigh);
    baselineH.time = tfaHigh.time;
    baselineH.dimord = tfaHigh.dimord;
    baselineH.powspctrm = repmat(baselineH.powspctrm, [1,1,length(baselineH.time), size(tfaHigh.powspctrm, 1)]);
    baselineH.powspctrm = permute(baselineH.powspctrm, [4,1,2,3]);
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    tfaHigh = ft_math(cfg, tfaHigh, baselineH);
    
    for freq=1:19
        for ch=1:length(d.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_h = squeeze(squeeze(tfaHigh.powspctrm(:,ch,freq,:)));
            betas_h(freq,ch,:,:) = design'\Y_h;
        end
    end
    
elseif strcmp(freqRange, 'low')
    % do the same for low frequency TFR
    cfg=[];
    cfg.latency = [-1+1/fs -0.25];
    cfg.avgoverrpt = 'yes';
    cfg.avgovertime = 'yes';
    baselineL = ft_selectdata(cfg, tfaLow);
    baselineL.time = tfaLow.time;
    baselineL.dimord = tfaLow.dimord;
    baselineL.powspctrm = repmat(baselineL.powspctrm, [1,1,length(baselineL.time), size(tfaLow.powspctrm, 1)]);
    baselineL.powspctrm = permute(baselineL.powspctrm, [4,1,2,3]);
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    tfaLow = ft_math(cfg, tfaLow, baselineL);
    
    for freq=1:15
        for ch=1:length(d.label);
            design = [ones(size(p1amp(ch,:))); p1amp(ch,:)];
            design(2,:) = (design(2,:)-mean(design(2,:)))./std(design(2,:));
            Y_l = squeeze(squeeze(tfaLow.powspctrm(:,ch,freq,:)));
            betas_l(freq,ch,:,:) = design'\Y_l;
        end
    end
end

%% planar gradiant transformation of beta weights
% make timelocked structure where planar gradient transformation can be
% applied to (that's why dimord is strange)
if strcmp(freqRange, 'high')
    tlh=[];
    tlh.avg    = squeeze(betas_h(:,:,2,:));
    tlh.time   = tfaHigh.time;
    tlh.dimord = 'subj_chan_time';
    tlh.label  = tfaHigh.label;
    tlh.grad   = tfaHigh.grad;
    %
    % % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tlh);
    cfg.planarmethod    = 'sincos';
    tlhPlanar           = ft_megplanar(cfg, tlh);
    
    % demeaning/baseline correction not possible on this data structure. Do
    % manual baseline correction. NOT NECESSARY. DOESN'T CHANGE ANYTHING
    %{
t3            = nearest(tlhPlanar.time, -0.6);
t4            = nearest(tlhPlanar.time, -0.1);
blcorr        = rmfield(tlhPlanar, 'trial');
blcorr.trial  = repmat(mean(tlhPlanar.trial(:,:,t3:t4), 3), [1,1,size(tlhPlanar.trial,3)]); % take average of baseline window
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'trial';
tlhPlanar = ft_math(cfg, tlhPlanar, blcorr);
    %}
    cfg           = [];
    bhPlanarCmb  = ft_combineplanar(cfg,tlhPlanar);
    
    % put it back in a freq-data structure
    bhPlanarCmb.powspctrm = permute(bhPlanarCmb.trial, [2,1,3]);
    bhPlanarCmb           = rmfield(bhPlanarCmb, 'trial');
    bhPlanarCmb.freq      = tfaHigh.freq;
    bhPlanarCmb.dimord    = 'chan_freq_time';
    
elseif strcmp(freqRange, 'low')
    % Do the same for low frequencies.
    tll=[];
    tll.avg    = betas_l;
    tll.time   = tfaLow.time;
    tll.dimord = 'subj_chan_time';
    tll.label  = tfaLow.label;
    tll.grad   = tfaLow.grad;
    %
    % % planar combination
    cfg                 = [];
    cfg.feedback        = 'no';
    cfg.method          = 'template';
    cfg.neighbours      = ft_prepare_neighbours(cfg, tll);
    cfg.planarmethod    = 'sincos';
    tllPlanar           = ft_megplanar(cfg, tll);
    
    cfg           = [];
    blPlanarCmb  = ft_combineplanar(cfg,tllPlanar);
    
    % put it back in a freq-data structure
    blPlanarCmb.powspctrm = permute(blPlanarCmb.trial, [2,1,3]);
    blPlanarCmb           = rmfield(blPlanarCmb, 'trial');
    blPlanarCmb.freq      = tfaLow.freq;
    blPlanarCmb.dimord    = 'chan_freq_time';
end
%% Normalize beta weights
% also doens't change much
if strcmp(freqRange, 'high')
    t5              = nearest(bhPlanarCmb.time, -0.6);
    t6              = nearest(bhPlanarCmb.time, -0.1);
    muH              = rmfield(bhPlanarCmb, 'powspctrm');
    muH.powspctrm    = nanmean(bhPlanarCmb.powspctrm(:,:,t5:t6), 3);
    muH.powspctrm    = repmat(muH.powspctrm, [1, 1, length(bhPlanarCmb.time)]);
    sigmaH           = rmfield(bhPlanarCmb, 'powspctrm');
    sigmaH.powspctrm = nanstd(bhPlanarCmb.powspctrm(:,:,t5:t6),[],3);
    sigmaH.powspctrm = repmat(sigmaH.powspctrm, [1,1, length(bhPlanarCmb.time)]);
    
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2)./x3';
    bhPlanarCmbZ = ft_math(cfg, bhPlanarCmb, muH, sigmaH);
    
elseif strcmp(freqRange, 'low')
    % % Do the same for low freq
    t7              = nearest(blPlanarCmb.time, -0.6);
    t8              = nearest(blPlanarCmb.time, -0.1);
    muL              = rmfield(blPlanarCmb, 'powspctrm');
    muL.powspctrm    = mean(blPlanarCmb.powspctrm(:,:,t7:t8), 3);
    muL.powspctrm    = repmat(muL.powspctrm, [1, 1, length(blPlanarCmb.time)]);
    sigmaL           = rmfield(blPlanarCmb, 'powspctrm');
    sigmaL.powspctrm = std(blPlanarCmb.powspctrm(:,:,t7:t8),[],3);
    sigmaL.powspctrm = repmat(sigmaL.powspctrm, [1,1, length(blPlanarCmb.time)]);
    
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2)./x3';
    blPlanarCmbZ = ft_math(cfg, blPlanarCmb, muL, sigmaL);
end
%% Save
if strcmp(freqRange, 'high')
    if isPilot
        filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_tfH_%d', subj, zeropoint);
    else
        filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tfH_%d', subj, zeropoint);
    end
    save(fullfile([filename '.mat']), 'betas_h','bhPlanarCmb', 'bhPlanarCmbZ','lat', '-v7.3');
elseif strcmp(freqRange, 'low')
    if isPilot
        filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_tfL_%d', subj, zeropoint);
    else
        filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_tfL_%d', subj, zeropoint);
    end
    save(fullfile([filename '.mat']), 'betas_l', 'blPlanarCmb','blPlanarCmbZ', 'lat', '-v7.3');
end

diary off
movefile(diaryname, fullfile([filename '.txt']));

