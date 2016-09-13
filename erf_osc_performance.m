function erf_osc_performance(subj, isPilot)

if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

datainfo_erf_osc; % load subject specific info.

cfg=[];
if isPilot
    cfg.dataset = pilotsubjects(subj).dataset;
    load(pilotsubjects(subj).logfile);% load log file
else
    cfg.dataset = subjects(subj).dataset;
    load(subjects(subj).logfile);
end
cfg.datafile = cfg.dataset;
cfg.headerfile = cfg.dataset;
cfg.trialfun = 'mytrialfun';
cfg.trialdef.prestim = min(log.realBaselineDuration, log.setBaselineDuration);
cfg.trialdef.poststim = log.completeDurationGrating;
cfg.catchtrial = log.trlNoShift;
cfg.continuous = 'yes';
cfg = ft_definetrial(cfg);


% preprocess data
cfg.channel = 'UPPT002';
data = ft_preprocessing(cfg);

sum(data2.trialinfo(:,5)==1)/600


end


