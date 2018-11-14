function erf_osc_analysis_erf_virtualchan(subj, erfoi)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(erfoi)
    erfoi = 'reversal';
end

ft_diary('on')
erf_osc_datainfo;


lcmvData = erf_osc_analysis_lcmv_orientation(subj, erfoi,'reverse', 'lp'); % optimal dipole orientation is chosen based on [0 0.5] after stimulus-reversal.
% data is timelocked to reversal or behavioral response, and low pass
% filtered. 
load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_gamma_virtual_channel.mat', subj, subj),'gammaPow'); % load gamma power

cfg=[];
cfg.vartrllength = 2;
if strcmp(erfoi,'reversal')
    cfg.preproc.baselinewindow = [-0.1 0];
elseif strcmp(erfoi, 'motor')
    cfg.preproc.baselinewindow = [-0.6 -0.5];
end
cfg.preproc.demean = 'yes';
cfg.keeptrials = 'yes';
tlck = ft_timelockanalysis(cfg, lcmvData);

[val idx] = sort(gammaPow, 'descend');
cfg=[];
if strcmp(erfoi,'reversal')
    cfg.latency = [0 0.5];
elseif strcmp(erfoi, 'motor')
    cfg.latency = [-0.5 0];
end
tlck = ft_selectdata(cfg, tlck);
qsize  = round(length(tlck.trialinfo)/4);
cfg1=[];
cfg1.trials = idx(1:qsize);
q1 = ft_selectdata(cfg, tlck);
cfg1.trials = idx(end-qsize+1:end);
q4 = ft_selectdata(cfg1, tlck);

tl_q1 = ft_timelockanalysis([], q1);
tl_q4 = ft_timelockanalysis([], q4);

if strcmp(erfoi,'motor')
    tlck=[];
end
%% save
filename = sprintf('/project/3011085.02/analysis/erf/sub-%03d/sub-%03d_erf_virtualchan_%s', subj, subj, erfoi);
save(fullfile([filename '.mat']), 'tlck', 'qsize', '-v7.3')
ft_diary('off')

end