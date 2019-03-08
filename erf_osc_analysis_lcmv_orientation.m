function [data] = erf_osc_analysis_lcmv_orientation(subj, erfoi, filtdir, dofilter, filtfreq)
% take the lcmv dat
if nargin<1 || isempty(subj)
    error('subject number was not specified')
end
if nargin<2 || isempty(erfoi)
    erfoi = 'reversal';
end
if nargin<3 || isempty(filtdir)
    filtdir = 'onepass-zerophase';
end
if nargin<4 || isempty(dofilter)
    dofilter = 'lp';
end
if nargin<5 || isempty(filtfreq)
    filtfreq = 30;
end

ft_diary('on')

erf_osc_datainfo;
load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel.mat',subj, subj), 'lcmvData');

cfg=[];
if strcmp(erfoi, 'reversal')
    cfg.offset = -(lcmvData.trialinfo(:,5)-lcmvData.trialinfo(:,4));
elseif strcmp(erfoi, 'motor')
    cfg.offset = -(lcmvData.trialinfo(:,6)-lcmvData.trialinfo(:,4));
end
data = ft_redefinetrial(cfg, lcmvData);

if strcmp(dofilter, 'bp')
    cfg=[];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = filtfreq;
    cfg.bpfiltdir = 'onepass-zerophase';
    cfg.bpfilttype = 'firws';
elseif strcmp(dofilter, 'lp')
    cfg=[];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = filtfreq;
    cfg.lpfilttype = 'firws';
    if strcmp(filtdir, 'reverse')
        cfg.lpfiltdir = 'onepass-reverse-zerophase';
    else
        cfg.lpfiltdir = filtdir;
    end
end
data = ft_preprocessing(cfg, data);

cfg=[];
cfg.vartrllength = 2;
if ~strcmp(erfoi, 'motor')
    cfg.preproc.baselinewindow = [-0.1 0];
end
cfg.preproc.demean = 'yes';
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, data);

cfg=[];
if strcmp(erfoi,'motor')
    cfg.latency = [-0.5 0];
else
    cfg.latency = [0 0.5];
end
tlck = ft_selectdata(cfg, tlck);

cfg=[];
cfg.avgoverrpt='yes';
tlck_avg = ft_selectdata(cfg, tlck);

[u,s,v]=svd(tlck_avg.trial, 'econ');

data.trial = u(:,1)'*data.trial; % take only first 'principle component' or
% for average of 2 components:dot(u(:,1)'*data.trial, u(:,2)'*data.trial,2)
data.label = {'gam_pow'};

ft_diary('off')
end