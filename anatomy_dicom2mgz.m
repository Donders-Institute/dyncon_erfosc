function anatomy_dicom2mgz(subj, cfg)
% anatomy_dicom2mgz takes the the subject info data structure (or subject string as 'sXX') 
%   
%   Picks up the dicom files, reslices the image and creates a .mgz file (spm coordsyst)

subjects = ft_getopt(cfg, 'subjects'); % if isPilot, cfg.subjects=pilotsubjects
mri = ft_read_mri(subjects(subj).mri); % read in the dicom files
mgz_dir = subjects(subj).mridir; % filename for saving
preproc_dir = fullfile(mgz_dir, '/preproc/');
if ~exist(preproc_dir)
    mkdir(preproc_dir)
end

% save the images in the mgz format
cfg             = [];
cfg.filename    = fullfile(mgz_dir, '/preproc/mri.mgz');
cfg.filetype    = 'mgz';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri);

end

