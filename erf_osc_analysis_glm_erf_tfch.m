function erf_osc_analysis_glm_erf_tfch(subj, isPilot, glm)
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
    glm = 'fixchan'; % other option is 'fixtime'
end
if isempty(glm)
    glm='fixchan';
end


%% Initiate Diary
workSpace = whos;
diaryname = tempname;
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

%%
% first six subjects, maximum channels
c(1)=match_str(data_dss.label,'MRO22');
c(2)=match_str(data_dss.label,'MRO52');
c(3)=match_str(data_dss.label,'MRO23');
c(4)=match_str(data_dss.label,'MRO52');
c(5)=match_str(data_dss.label,'MRO11');
c(6)=match_str(data_dss.label,'MRO32');
lat = [0.09917, 0.06333, 0.1092, 0.07417, 0.07, 0.11];

if strcmp(glm, 'fixtime')
    %% all channels, one time
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.baselinetype = 'absolute';
    cfg.baseline = [-0.6 -0.1];
    h = ft_freqbaseline(cfg, tfaHigh);
    l = ft_freqbaseline(cfg, tfaLow);
    
    cfg=[];
    cfg.latency=[0.05 0.15];
    cfg.avgovertime='yes';
    l = ft_selectdata(cfg, l);
    h = ft_selectdata(cfg, h);
    
    cfg=[];
    cfg.channel = data_dss.label(c(subj));
    d  =ft_selectdata(cfg, data_dss);
    d.trial = cat(1,d.trial{:});
    
    idx = nearest(d.time{1}, lat(subj));
    p1 = d.trial(:,idx)'; % p1 amplitude per trial
    
    design = [ones(size(p1)); p1];
    
    for freq=1:19
        Y_h = (squeeze(h.powspctrm(:,:,freq)));
        betas_h(:,:,freq) = design'\Y_h;
    end
    bh = h;
    bh.powspctrm = betas_h(2,:,:); % save betas in powspctrm structure, in right format for later use of ft_freqdescriptives.
    bh.dimord='rpt_chan_freq';
    
    % do the same for low frequency TFR
    for freq=1:15
        Y_l = squeeze(l.powspctrm(:,:,freq));
        betas_l(:,:,freq) = design'\Y_l;
    end
    bl = l;
    bl.powspctrm = betas_l(2,:,:);
    bl.dimord='rpt_chan_freq';
    
    
elseif strcmp(glm, 'fixchan')
    %% one channel, over time
    select data only for maximum channel
    cfg=[];
    cfg.channel = data_dss.label(c(subj));
    l = ft_selectdata(cfg, tfaLow);
    h = ft_selectdata(cfg, tfaHigh);
    d  =ft_selectdata(cfg, data_dss);
    d.trial = cat(1,d.trial{:});
    
    idx = nearest(d.time{1}, lat(subj));
    p1 = d.trial(:,idx)'; % p1 amplitude per trial
    
    design = [ones(size(p1)); p1];
    
    baseline correct TFR
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.baselinetype = 'absolute';
    cfg.baseline = [-0.6 -0.1];
    h = ft_freqbaseline(cfg, h);
    l = ft_freqbaseline(cfg, l);
    
    for freq=1:19
        Y_h = squeeze(squeeze(h.powspctrm(:,:,freq,:)));
        betas_h(:,1,:,freq) = design'\Y_h;
    end
    bh = h;
    bh.powspctrm = permute(betas_h(2,:,:,:), [1,2,4,3]); % save betas in powspctrm structure, in right format for later use of ft_freqdescriptives.
    bh.label = {'chan01'};
    
    % do the same for low frequency TFR
    for freq=1:15
        Y_l = squeeze(squeeze(l.powspctrm(:,:,freq,:)));
        betas_l(:,1,:,freq) = design'\Y_l;
    end
    bl = l;
    bl.powspctrm = permute(betas_l(2,:,:,:), [1,2,4,3]);
    bl.label={'chan01'};
    
    %% Save
    if isPilot
        if strcmp(glm, 'fixtime')
            filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_erf_fch', subj);
        elseif strcmp(glm, 'fixchan')
            filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/glm_erf_tf', subj);
        end
    else
        if strcmp(glm, 'fixtime')
            filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_erf_fch', subj);
        elseif strcmp(glm, 'fixch')
            filename = sprintf('/project/3011085.02/results/erf/sub-%03d/glm_erf_tf', subj);
        end
    end
    save(fullfile([filename '.mat']),'betas_h', 'betas_l', 'bl', 'bh','-v7.3');
    diary off
    movefile(diaryname, fullfile([filename '.txt']));
    
