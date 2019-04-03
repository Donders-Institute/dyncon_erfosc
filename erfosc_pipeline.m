%% Set up
% setting up the environment by running the startup script.
startup

erf_osc_datainfo; % get subject specific information

dosave = false;

%% preprocessing
% seperately preprocess the MRI data using fieldtrip, freesurfer and
% workbench. Do coregistration, segmentation and create 2D and 3D
% sourcemodels and a headmodel. This script cannot be automatized and
% should be done interactively by walking through the steps.
for subj=allsubs
    erfosc_execute_pipeline('erf_osc_preprocessing_mri', subj);
end



% MEG preprocessing is done semi-automatically; bad channels/trials,
% jump/muscle and eye blink artifacts have to be manually verified.
% ICA decomposition should be inspected and components to be deleted should
% be specified in erf_osc_datainfo, under subjects(subj).ecgcomp/eyecomp
% before continuing. This function saves the multiple files: one containing
% sample information on artifacts, one with ICA spatial filters, one with
% headmotion data and one with the preprocessed MEG data.
    global ft_default;
    ft_default = [];
    ft_default.reproducescript = datestr(now, 30);
for subj=allsubs
    erf_osc_preprocessing_artifact(subj, false, true, 5, dosave)

% seperately, preprocess the eyetracker data. Transform X- and Y-gaze
% positions to visual degrees relative to fixation. Save it on disk
% together with pupil diameter data, both for data time locked to stimulus
% onset and stimulus change.

    [eyedata_onset, eyedata_shift] = erfosc_preprocessing_eyedata(subj, dosave);

    %% Single Subject Analysis
% Reaction times
    erf_osc_analysis_rt(subj, dosave);

% Stimulus induced effects

    erf_osc_analysis_tfa(subj, false, 'low', 'onset', dosave) % low frequencies
    erf_osc_analysis_tfa(subj, false, 'high', 'onset', dosave) % high frequencies

%% power spectra
% compute channel level induced power spectra in low and high frequencies.
    erf_osc_analysis_pow(subj, false);

%% 2D source parcels    
    erfosc_execute_pipeline('erfosc_script', subj, {'dolcmv_parc', true}, {'dolcmv_norm', true}, {'savefreq', dosave});
    
%% Source power contrast
    erfosc_execute_pipeline('erfosc_script', subj, {'dodics_gamma', true}, {'dodics_lowfreq', true}, {'savefreq', dosave});

%% correlation of rt with gamma in 3D sourcespace
    erfosc_execute_pipeline('erfosc_script', subj, {'docorr_gamma_rt', true}, {'savefreq', dosave});

%% Correlation Gamma ERF
    erfosc_execute_pipeline('erfosc_script', subj, {'docorrpow_lcmv',1}, {'dofreq_short', 1}, {'savefreq', dosave}, {'dolcmv_norm', 1});
    erfosc_execute_pipeline('erfosc_script', subj, {'docorrpow_lcmv_lowfreq',1}, {'dofreq_short_lowfreq', 1}, {'savefreq', dosave});
    ft_default.reproducescript = [];
end

%% Group Analysis
%% Induced power
erf_osc_analysis_tfa_allsubs('low', 'onset', dosave);
erf_osc_analysis_tfa_allsubs('high', 'onset', dosave);

%% Correlation with RT
% computing correlations and partial correlations between jitter, gamma
% power and rt.
erf_osc_analysis_corr([],[],'gamma_rt')

%% Correlation Gamma ERF
erfosc_execute_pipeline('erfosc_script', 1, {'dostat_pow_erf',1}, {'GA', 1}, {'whichFreq', 1}, {'savefreq', dosave});
erfosc_execute_pipeline('erfosc_script', 1, {'dostat_pow_erf',1}, {'GA', 1}, {'whichFreq', 2},, {'savefreq', dosave});


%% Correlation ERF RT
erfosc_execute_pipeline('erfosc_script', 1, {dostat_erf_rt, 'true'});

%% Plotting
plot_figures_chapter

ft_default.reproducescript = [];
