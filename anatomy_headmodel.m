function anatomy_headmodel(subj, cfg)

subjects              = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir   = subjects(subj).mridir;
anatomy_preproc_dir   = fullfile(anatomy_preproc_dir, '/preproc/');
headmodel_filename    = fullfile(anatomy_preproc_dir, 'headmodel.mat');

mni_resliced_filename       = fullfile(anatomy_preproc_dir, 'mni_resliced.mgz');
transform                   = fullfile(anatomy_preproc_dir, 'transform_vox2ctf.mat');
load(transform);

mri                         = ft_read_mri(mni_resliced_filename);
mri.coordsys                = 'ctf';
mri.transform               = transform_vox2ctf;

cfg = [];
cfg.output = 'brain';
seg = ft_volumesegment(cfg, mri);

cfg = [];
cfg.method = 'projectmesh';
cfg.numvertices = 10000;
bnd = ft_prepare_mesh(cfg, seg);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, bnd);

save(headmodel_filename, 'headmodel');

end

