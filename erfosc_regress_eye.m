function residual = erfosc_regress_eye(eyedata, chan, confound)
% regresses out confound from the data. Optimized for eye tracking data.
% Residuals of the regressionanalysis are given as output
%
% INPUT
%   eyedata: fieldtrip data structure
%   chan: from which channels the confound should be regressed out
%   confound: the channel names of the confounding variables
%
% OUTPUT
%   residual: residual after ft_regressconfound

cfg=[];
cfg.channel = chan;
cfg.vartrllength=2;
cfg.keeptrials = 'yes';
pupilpre = ft_timelockanalysis(cfg, eyedata);

cfg.channel = confound;
cc = ft_timelockanalysis(cfg, eyedata);

cfg=[];
cfg.output = 'residual';
for k=1:numel(confound)
    idx = match_str(cc.label, confound{k});
    tmp(:,k) = nanmean(squeeze(cc.trial(:,k,:)),2);
end
cfg.confound = tmp;
residual = ft_regressconfound(cfg, pupilpre);
