function erf_osc_analysis_glm_erf_tfch(subj, isPilot)
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
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/tfa.mat', subj));
    load(sprintf('/project/3011085.02/results/erf/sub-%03d/dss.mat', subj), 'data_dss');
    %     load(sprintf('/project/3011085.02/results/erf/sub-%03d/timelock.mat', subj));
end
fs=data_dss.fsample;
nTrials = length(data_dss.trial);

%% select p1 window for regression

tlck = ft_timelockanalysis([], data_dss);
t1 = nearest(tlck.time, 0.06);
t2 = nearest(tlck.time, 0.12);
cfg=[];
cfg.latency = [tlck.time(t1) tlck.time(t2)];
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

% first six subjects, maximum channels
% c(1)=match_str(data_dss.label,'MRO22');
% c(2)=match_str(data_dss.label,'MRO52');
% c(3)=match_str(data_dss.label,'MRO23');
% c(4)=match_str(data_dss.label,'MRO52');
% c(5)=match_str(data_dss.label,'MRO11');
% c(6)=match_str(data_dss.label,'MRO32');
% lat = [0.09917, 0.06333, 0.1092, 0.07417, 0.07, 0.11];

 %% Regression p1 amplitude over time-frequency-channel
    %select data only for maximum channel
    cfg=[];
    cfg.channel = 'all';%data_dss.label(c(subj));
    l = ft_selectdata(cfg, tfaLow);
    h = ft_selectdata(cfg, tfaHigh);
    d  =ft_selectdata(cfg, data_dss);
    d.trial = cat(3,d.trial{:});
    
    % baseline correct TFR
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.baselinetype = 'relchange';
    cfg.baseline = [-0.6 -0.1];
    h = ft_freqbaseline(cfg, h);
    l = ft_freqbaseline(cfg, l);
    
    idx = [nearest(d.time{1}, lat(1)), nearest(d.time{1}, lat(2))];
    p1 = squeeze(mean(d.trial(:,idx,:),2));
    
    for freq=1:19
        for ch=1:length(d.label);
            design = p1(ch,:);
            Y_h = squeeze(squeeze(h.powspctrm(:,ch,freq,:)));
            betas_h(freq,ch,:) = design'\Y_h;
        end
    end
    bh=h;
    bh.powspctrm = permute(betas_h, [2,1,3]);
    bh.label=h.label;
    bh.dimord = 'chan_freq_time';
    
    % do the same for low frequency TFR
    for freq=1:15
        for ch=1:length(d.label);
            design = p1(ch,:);
            Y_l = squeeze(squeeze(l.powspctrm(:,ch,freq,:)));
            betas_l(freq,ch,:) = design'\Y_l;
        end
        bl = l;
        bl.powspctrm = permute(betas_l, [2,1,3]);
        bl.label=l.label;
        bl.dimord = 'chan_freq_time';
    end
    
    %% planar gradiant transformation of beta weights
tl=[];
tl.avg    = squeeze(betas_trials(2,:,:))';
tl.time   = data_dss.time{1};
tl.dimord = 'chan_time';
tl.label  = data_dss.label;
tl.grad   = data_dss.grad;
% 
% % planar combination
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.neighbours      = ft_prepare_neighbours(cfg, tl);
cfg.planarmethod    = 'sincos';
tlPlanar            = ft_megplanar(cfg, tl);
cfg                 = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.25 0];
tlPlanarCmb         = ft_combineplanar(cfg,tlPlanar);
%     

%%
t1        = nearest(tl.time, -2);
t2        = nearest(tl.time, 0);
mu        = rmfield(tl, 'avg');
mu.avg    = mean(tl.avg(:,t1:t2), 2);
mu.avg    = repmat(mu.avg, [1, length(tl.time)]);
sigma     = rmfield(tl, 'avg');
sigma.avg = std(tl.avg(:,t1:t2),[],2);
sigma.avg = repmat(sigma.avg, [1, length(tl.time)]);

cfg=[];
cfg.parameter = 'avg';
cfg.operation = '(x1-x2)./x3';
tlPlanarCmbZ2 = ft_math(cfg, tl, mu, sigma);

    %% Save
    if isPilot
            filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_erf_tf', subj);
    else
            filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_erf_tf', subj);
    end
    save(fullfile([filename '.mat']),'betas_h', 'betas_l', 'bl', 'bh', 'betas_h','lat', '-v7.3');
    diary off
    movefile(diaryname, fullfile([filename '.txt']));
    
