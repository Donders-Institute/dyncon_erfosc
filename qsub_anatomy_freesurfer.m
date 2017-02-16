function qsub_anatomy_freesurfer(cfg, subj, script_number)



if script_number==1
    % freesurfer script 1
    shell_script = '/project/3011085.02/scripts/erfosc/anatomy_freesurfer.sh';
elseif script_number==2
    shell_script = '/project/3011085.02/scripts/erfosc/anatomy_freesurfer2.sh';
elseif script_number==3
    shell_script = '/project/3011085.02/scripts/erfosc/anatomy_postfreesurferscript.sh';
end

mri_dir = cfg.subjects(subj).mridir;
preproc_dir = 'preproc';

% create the string that is executed in the linux terminal
command = [shell_script, ' ', mri_dir, ' ', preproc_dir];

% call the script
system(command);


