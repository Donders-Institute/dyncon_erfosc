function anatomy_mgz2ctf(subj, cfg)
%streams_anatomy_mgz2ctf reads in the resliced volume created with
%streams_anatomy_mgz2mni and the transformation matrix and calls
%ft_volumerealign to perform interactive coregistration to CTF coordinate
%headspace. It then saves the transformation matrix


%% Initialize the variables
subjects                  = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir       = subjects(subj).mridir;
anatomy_preproc_dir       = fullfile(anatomy_preproc_dir, '/preproc/');
resliced_filename         = fullfile(anatomy_preproc_dir, 'mni_resliced.mgz');
transformation_matrix_mni = fullfile(anatomy_preproc_dir, 'transform_vox2mni');

%read in the resliced volume
mri_resliced_mni = ft_read_mri(resliced_filename);

% check if the transformation matrix exist, if so read it in 
if exist(transformation_matrix_mni, 'file');

  load(transformation_matrix_mni);
  
  if isequal(mri_resliced_mni.transform, transformation_matrix_mni)
    % do nothing
  else
    fprintf('adding coregistration information to the mrifile of subject %s\n', subj);

    mri_resliced_mni.transform = transform_vox2mni;

    cfg = [];
    cfg.parameter = 'anatomy';
    cfg.filename  = resliced_filename;
    cfg.filetype  = 'mgz';
    ft_volumewrite(cfg, mri_resliced_mni);

    clear mri;
  end

end
%% Interactive realignment to the CTF conventions
% do an initial coregistration
cfg             = [];
cfg.interactive = 'yes';
cfg.coordsys    = 'ctf';
mri = ft_volumerealign(cfg, mri_resliced_mni);

% refine the coregistration based on the polhemus point cloud
polhemus = ft_read_headshape(subjects(subj).headshape);
polhemus.unit='cm';

cfg                         = [];
cfg.coordsys                = 'ctf';
cfg.parameter               = 'anatomy';
cfg.viewresult              = 'yes';
cfg.method                  = 'headshape';
cfg.headshape.headshape     = polhemus;
cfg.headshape.interactive   = 'yes';
cfg.headshape.icp           = 'yes';
mri_ctf    = ft_volumerealign(cfg, mri_resliced_mni);


% save the transformation matrix
transform_vox2ctf = mri_ctf.transform;
filename_vox2ctf  = fullfile(anatomy_preproc_dir, 'transform_vox2ctf');
save(filename_vox2ctf, 'transform_vox2ctf');
  
end

