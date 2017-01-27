
% load structural MRI and polhemus headshape. Realign the structural scan
% to CTF coordinate system and segment.

subj = 1;
isPilot = false;
preproc = '2d';



if strcmp(preproc, '2d') % cortical sheet
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CORTICAL SHEET, FREESURFER, PREPROCESSING %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    erf_osc_datainfo;
    if isPilot
        cfg.subjects = pilotsubjects;
    else
        cfg.subjects = subjects;
    end
    
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
        'memreq', 1024^3 * 6, 'timreq', 720*60, 'batchid', 'erfosc_freesurfer1')
    
    %% Check-up and white matter segmentation cleaning if needed
    
    anatomy_volumetricQC(subj, cfg)
    
    anatomy_wmclean(subj, cfg)
    
    %% Freesurfer qsub2
    if ~ft_hastoolbox('qsub',1)
        addpath /project/3011085.02/scripts/fieldtrip/qsub;
    end
    
    qsubfeval(@qsub_anatomy_freesurfer, cfg, subj, 2,...
        'memreq', 1024^3 * 7, 'timreq', 720*60, 'batchid', 'erfosc_freesurfer2');
    
    %% Post-processing Freesurfer script: workbench HCP tool
    if ~ft_hastoolbox('qsub',1)
        addpath /project/3011085.02/scripts/fieldtrip/qsub;
    end
    
    qsubfeval(@qsub_anatomy_freesurfer, cfg, subj, 3,...
        'memreq', 1024^3 * 6, 'timreq', 480*60, 'batchid', 'erfosc_workbench');
    
    % Sourcemodel
    anatomy_sourcemodel2d(subject, cfg);
    
    % Headmodel
    anatomy_headmodel(subject, cfg);
    
    % Coregistration check
    anatomy_coregistration_qc(subject, cfg);
    
    
    
    
    
    
elseif strcmp(preproc, '3d') % normal, volumetric
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STANDARD (VOLUMETRIC) MRI PREPROCESSING %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% initiate diary
    workSpace = whos;
    diary('tmpDiary') % save command window output
    fname = mfilename('fullpath')
    datetime
    
    fid = fopen(fullfile([fname '.m']));
    tline = fgets(fid); % returns first line of fid
    while ischar(tline) % at the end of the script tline=-1
        disp(tline) % display tline
        tline = fgets(fid); % returns the next line of fid
    end
    fclose(fid);
    
    for i = 1:numel(workSpace) % list all workspace variables
        workSpace(i).name % list the variable name
        printstruct(eval(workSpace(i).name)) % show its value(s)
    end
    
    %% Load data, realign and segment
    erf_osc_datainfo;
    
    mri = ft_read_mri(pilotsubjects(subj).mri);
    mri.coordsys = 'mni';
    polhemus = ft_read_headshape(pilotsubjects(subj).headshape);
    polhemus.unit='cm';
    
    % align mri with MEG data through polhemus headshape
    cfg                         = [];
    cfg.parameter               = 'anatomy';
    cfg.viewresult              = 'yes';
    cfg.method                  = 'headshape';
    cfg.headshape.headshape     = polhemus;
    cfg.headshape.interactive   = 'yes';
    cfg.headshape.icp           = 'yes';
    cfg.coordsys                = 'ctf';
    mri2                         = ft_volumerealign(cfg, mri);
    
    % segment mri
    cfg             = [];
    cfg.write       = 'no';
    cfg.viewresult  = 'yes';
    [segmentedmri]  = ft_volumesegment(cfg, mri);
    
    %% Save
    filename = pilotsubjects(subj).segmentedmri;
    save(fullfile([filename '.mat']),'segmentedmri')
    diary off
    movefile('tmpDiary', fullfile([filename, '.txt']));
end

