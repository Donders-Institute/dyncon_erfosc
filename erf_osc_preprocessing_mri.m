function erf_osc_preprocessing_mri(subj, isPilot)
% load structural MRI and polhemus headshape. Realign the structural scan
% to CTF coordinate system and segment.


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

% align mri with MEG data through polhemus headshape
cfg                         = [];
cfg.parameter               = 'anatomy';
cfg.viewresult              = 'yes';
cfg.method                  = 'headshape';
cfg.headshape.headshape     = polhemus;
cfg.headshape.interactive   = 'yes';
cfg.headshape.icp           = 'yes';
cfg.coordsys                = 'ctf';
mri                         = ft_volumerealign(cfg, mri);

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