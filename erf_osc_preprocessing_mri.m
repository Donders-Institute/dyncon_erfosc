
% load structural MRI and polhemus headshape. Realign the structural scan
% to CTF coordinate system and segment.

if ~exist('subj'); subj = 1; end
if ~exist('isPilot'); isPilot = false; end

erf_osc_datainfo;
if isPilot
    cfg.subjects = pilotsubjects;
else
    cfg.subjects = subjects;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORTICAL SHEET, FREESURFER, PREPROCESSING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% converting dicoms to mgz format
anatomy_dicom2mgz(subj, cfg);

% reslicing to freesufer-friendly 256x256x256
anatomy_mgz2mni(subj, cfg);

anatomy_mgz2ctf(subj, cfg);

% Skullstriping
anatomy_skullstrip(subj, cfg);

%% Freesurfer scripts (creates subject-specific subdirectory in the directory where previous files are stored)
if ~ft_hastoolbox('qsub',1)
    addpath /project/3011085.02/scripts/fieldtrip/qsub;
end

qsubfeval(@qsub_anatomy_freesurfer, cfg, subj, 1,...
    'memreq', 1024^3 * 6, 'timreq', 720*60, 'batchid', sprintf('erfosc_freesurfer1_%d', subj))

%% Check-up and white matter segmentation cleaning if needed

anatomy_volumetricQC(subj, cfg)

anatomy_wmclean(subj, cfg)

%% Freesurfer qsub2
if ~ft_hastoolbox('qsub',1)
    addpath /project/3011085.02/scripts/fieldtrip/qsub;
end

qsubfeval(@qsub_anatomy_freesurfer, cfg, subj, 2,...
    'memreq', 1024^3 * 7, 'timreq', 720*60, 'batchid', sprintf('erfosc_freesurfer2_%d', subj));

%% Post-processing Freesurfer script: workbench HCP tool
if ~ft_hastoolbox('qsub',1)
    addpath /project/3011085.02/scripts/fieldtrip/qsub;
end

qsubfeval(@qsub_anatomy_freesurfer, cfg, subj, 3,...
    'memreq', 1024^3 * 6, 'timreq', 480*60, 'batchid', sprintf('erfosc_workbench_%d', subj));
%%
% Sourcemodel
anatomy_sourcemodel2d(subj, cfg);
anatomy_sourcemodel3d(subj, cfg); % specify resolution in cfg.resolution (in mm)

% Headmodel
anatomy_headmodel(subj, cfg);

% Coregistration check
anatomy_coregistration_qc(subj, cfg);




