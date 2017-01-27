function anatomy_workbench(subj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Strings for the command
shell_script      = '/home/language/kriarm/pro/streams/code/streams/matlab/streams_anatomy_postfreesurferscript.sh';
mri_dir           = '/home/language/kriarm/pro/streams/data/MRI/preproc';
subject_dir       = subj;

% streams_anatomy_freesurfer2.sh
command = [shell_script, ' ', mri_dir, ' ', subject_dir];

system(command);

end

