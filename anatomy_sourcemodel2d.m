function [sourcemodel] = anatomy_sourcemodel2d(subj, cfg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

subjects              = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir   = subjects(subj).mridir;
anatomy_preproc_dir   = fullfile(anatomy_preproc_dir, '/preproc/');
workbench_dir         = fullfile(anatomy_preproc_dir, '/workbench/');
sourcemodel2d_filename= fullfile(anatomy_preproc_dir, 'sourcemodel2d.mat'); %string for saving the sourcemodel file


% load in the cortical sheet
filename = fullfile(workbench_dir, ['preproc', '.L.midthickness.8k_fs_LR.surf.gii']);
filename2 = strrep(filename, '.L.', '.R.');

sourcemodel = ft_read_headshape({filename, filename2});

% get the necessary coregistration information
datapath = fullfile(anatomy_preproc_dir);
load(fullfile(datapath, 'transform_vox2mni'));
T1 = transform_vox2mni;
load(fullfile(datapath, 'transform_vox2ctf'));
T2 = transform_vox2ctf;

sourcemodel = ft_transform_geometry((T2/T1), sourcemodel);
sourcemodel.inside = sourcemodel.atlasroi>0;
sourcemodel = rmfield(sourcemodel, 'atlasroi');

save(sourcemodel2d_filename, 'sourcemodel');

end

