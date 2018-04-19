function [sourcemodel] = anatomy_sourcemodel3d(subj, cfg)

% MOUS_ANATOMY_SOURCEMODEL3D computes a 3D regular grid with 
% specified resolution based on an inverse warp of a template
% grid in MNI space.
%
% Changelog: 26-02-2013: use templates in MOUS/meg/templates/ directory.
% This means recomputing all sourcemodels for all subjects. Reason: the
% default sourcemodels in FieldTrip are not fully covering the top of the
% brain.
resolution = ft_getopt(cfg, 'resolution', 6);
fname   = fullfile('/project/3011085.02/scripts/fieldtrip/', 'template', 'sourcemodel', ['standard_sourcemodel3d',num2str(resolution),'mm.mat']);
load(fname);

subjects                  = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir       = subjects(subj).mridir;
anatomy_preproc_dir       = fullfile(anatomy_preproc_dir, '/preproc/');
resliced_filename         = fullfile(anatomy_preproc_dir, 'mni_resliced.mgz');
transformation_matrix_ctf = fullfile(anatomy_preproc_dir, 'transform_vox2ctf');
sourcemodel3d_filename      = fullfile(anatomy_preproc_dir, 'sourcemodel3d.mat');

% read in the resliced volume
mri = ft_read_mri(resliced_filename);
load(transformation_matrix_ctf);
mri.coordsys                = 'ctf';
mri.transform               = transform_vox2ctf;

% segment mri
% cfg = [];
% cfg.output = 'brain';
% segmentedmri = ft_volumesegment(cfg, mri);

% create the grid
cfg = [];
cfg.grid.warpmni    = 'yes';
cfg.grid.template   = sourcemodel;
cfg.grid.nonlinear  = 'yes';
cfg.mri = mri; %segmentedmri;
sourcemodel = ft_prepare_sourcemodel(cfg);
sourcemodel = ft_convert_units(sourcemodel, 'mm');

% remove the mri-structure from grid.cfg
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'mri');
sourcemodel.cfg = rmfield(sourcemodel.cfg, 'callinfo');

save(sourcemodel3d_filename, 'sourcemodel');
end