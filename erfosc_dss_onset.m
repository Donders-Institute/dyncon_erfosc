function [comp_onset, comp_shift, comp_shift2] = erfosc_dss_onset(data_onset, data_shift)

fs = data_onset.fsample;

% compute a spatial filter based on the onset-related ERF
cfg         = [];
cfg.latency = [-0.1 0.5-1/fs];
tmp_onset   = ft_selectdata(cfg, data_onset);
tmp_shift   = ft_selectdata(cfg, data_shift);

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq   = 60;
cfg.lpfilttype = 'firws';
tmp_onset = ft_preprocessing(cfg, tmp_onset);
tmp_shift = ft_preprocessing(cfg, tmp_shift);

params.time   = tmp_onset.time;
params.demean = 'prezero';

cfg                   = [];
cfg.method            = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params   = params;
cfg.dss.wdim          = 100;
cfg.numcomponent      = 10;
cfg.cellmode          = 'yes';
comp_onset            = ft_componentanalysis(cfg, tmp_onset);
params.time           = tmp_shift.time;
comp_shift            = ft_componentanalysis(cfg, tmp_shift);
%clear tmp_onset tmp_shift;

cfg           = [];
%cfg.topo      = comp_onset.topo;
cfg.unmixing  = comp_onset.unmixing;
cfg.topolabel = comp_onset.topolabel;
%comp_onset    = ft_componentanalysis(cfg, data_onset);
%comp_shift2   = ft_componentanalysis(cfg, data_shift);
comp_shift2    = ft_componentanalysis(cfg, tmp_shift);

%cfg.topo      = comp_shift.topo;
%cfg.unmixing  = comp_shift.unmixing;
%cfg.topolabel = comp_shift.topolabel;
%comp_shift    = ft_componentanalysis(cfg, data_shift);
