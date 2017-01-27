function anatomy_skullstrip(subj, cfg)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Initialize the variables

subjects                  = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir       = subjects(subj).mridir;
anatomy_preproc_dir       = fullfile(anatomy_preproc_dir, '/preproc/');
resliced_mni_filename     = fullfile(anatomy_preproc_dir, 'mni_resliced.mgz');
mri_skullstrip            = fullfile(anatomy_preproc_dir, 'skullstrip');

% read in the .mgz file created with anatomy_mgz2mni
mri_resliced_mni             = ft_read_mri(resliced_mni_filename);

% FSL variables
threshold       = 0.5;
T               = inv(mri_resliced_mni.transform);
center          = round(T(1:3,4))';

% name for the temporary nifti file
temp   = fullfile(anatomy_preproc_dir, 'nifti_tmp');

%% Skullstrip via FSL

% Convert to nifti temporarily and save;
cfg = [];
cfg.filename = temp;
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_resliced_mni);

% Create the FSL command-string
str = ['/opt/fsl/5.0.9/bin/bet ',temp,'.nii ',temp];
str = [str,'-R -f ',num2str(threshold),' -c ', num2str(center),' -g 0 -m -v'];

% Call the FSL command-string
system(str);

% Read the FSL-based segmentation
seg  = ft_read_mri([temp,'-R.nii.gz']);
delete([temp,'.nii']);
delete([temp,'-R.nii.gz']);
delete([temp,'-R_mask.nii.gz']);

% Save the FSL-based segmentation in .mgz
cfg = [];
cfg.filename = mri_skullstrip;
cfg.filetype = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, seg);

% Check the plot already now
skullstrip = ft_read_mri([mri_skullstrip '.mgz']);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, skullstrip);

end

