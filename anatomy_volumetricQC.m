function anatomy_volumetricQC(subj, cfg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

subjects                  = ft_getopt(cfg, 'subjects');
anatomy_preproc_dir       = subjects(subj).mridir;
anatomy_preproc_dir       = fullfile(anatomy_preproc_dir, '/preproc/');

t1              = fullfile(anatomy_preproc_dir, 'mri', 'T1.mgz'); % 'mri' needed?
normalization2  = fullfile(anatomy_preproc_dir, 'mri', 'brain.mgz');
white_matter    = fullfile(anatomy_preproc_dir, 'mri', 'wm.mgz');
white_matter_old = fullfile(anatomy_preproc_dir, 'mri', 'wm_old.mgz');

% Show T1
mri = ft_read_mri(t1);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subj ' ' 'T1'], 'numbertitle', 'off');

% Show skullstripped image
mri = ft_read_mri(normalization2);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subj ' ' 'skull-stripped'], 'numbertitle', 'off');

% Show white matter image
mri = ft_read_mri(white_matter);
cfg = [];
cfg.interactive = 'yes';
ft_sourceplot(cfg, mri);
set(gcf, 'name', [subj ' ' 'white matter'], 'numbertitle', 'off');

if exist(white_matter_old)
  
  mri = ft_read_mri(white_matter_old);
  cfg = [];
  cfg.interactive = 'yes';
  ft_sourceplot(cfg, mri);
  set(gcf, 'name', [subj ' ' 'white matter old'], 'numbertitle', 'off');
  
end

end

