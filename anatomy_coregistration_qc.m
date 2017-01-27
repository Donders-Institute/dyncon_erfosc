function anatomy_coregistration_qc(subj, cfg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subjects              = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir   = subjects(subj).mridir;
anatomy_preproc_dir   = fullfile(anatomy_preproc_dir, '/preproc/');
headmodel   = fullfile(anatomy_preproc_dir, 'headmodel.mat');
sourcemodel = fullfile(anatomy_preproc_dir, 'sourcemodel.mat');

load(headmodel)
load(sourcemodel)

figure; hold on;
ft_plot_vol(headmodel, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;

end

