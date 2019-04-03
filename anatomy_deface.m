function anatomy_deface(subj, cfg)

ctfspace = ft_getopt(cfg, 'ctfspace');
subjects               = ft_getopt(cfg, 'subjects');
anatomy_dir    = subjects(subj).mridir;

if ctfspace
    preproc_dir            = fullfile(anatomy_dir, '/preproc/');
    mni_resliced_filename  = fullfile(preproc_dir, 'mni_resliced.mgz');
    filename_vox2ctf       = fullfile(preproc_dir, 'transform_vox2ctf.mat');
    s = 'mri_ctf';
    
    
    load(filename_vox2ctf)
    mri = ft_read_mri(mni_resliced_filename);
    mri.transform = eye(4);
    mri = ft_transform_geometry(transform_vox2ctf, mri);
    
else
    s = 'mri';
    mri = ft_read_mri(subjects(subj).mri);
end
cfg=[];
cfg.method = 'spm';
mri = ft_defacevolume(cfg, mri);
ft_defacevolume([], mri);

save([anatomy_dir, sprintf('sub-%03d_%s_defaced.mat', subj, s)], 'mri');