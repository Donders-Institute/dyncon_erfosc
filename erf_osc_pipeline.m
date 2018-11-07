%% Set up
% setting up the environment by running the startup script.
startup

erf_osc_datainfo; % get subject specific information

%% preprocessing
% preprocessing is done semi-automatically; bad channels/trials,
% jump/muscle and eye blink artifacts have to be manually verified.
% ICA decomposition should be inspected and components to be deleted should
% be specified in erf_osc_datainfo, under subjects(subj).ecgcomp/eyecomp
% before continuing. This function saves the multiple files: one containing
% sample information on artifacts, one with ICA spatial filters, one with
% headmotion data and one with the preprocessed MEG data.
for subj=allsubs
    erf_osc_preprocessing_artifact(subj, false, false, 5)
end

% seperately, preprocess the eyetracker data. Transform X- and Y-gaze
% positions to visual degrees relative to fixation. Save it on disk
% together with pupil diameter data, both for data time locked to stimulus
% onset and stimulus change.
for subj=allsubs
    erf_osc_preprocessing_eyedata(subj);
end

% seperately preprocess the MRI data using fieldtrip, freesurfer and
% workbench. Do coregistration, segmentation and create 2D and 3D
% sourcemodels and a headmodel. This script cannot be automatized and
% should be done interactively by walking through the steps.
for subj=allsubs
    erfosc_execute_pipeline('erf_osc_preprocessing_mri', subj);
end



%% Analysis
%% Reaction times
for subj=allsubs
    erf_osc_analysis_rt(subj);
end


%% Stimulus induced effects
for subj=allsubs
    erf_osc_analysis_tfa(subj, false, 'low', 'onset') % low frequencies
    erf_osc_analysis_tfa(subj, false, 'high', 'onset') % high frequencies
end

% combine them for all subjects and make a relative contrast with the
% average baseline
erf_osc_analysis_tfa_allsubs('low', 'onset');
erf_osc_analysis_tfa_allsubs('high', 'onset');


%% power spectra
% compute channel level induced power spectra in low and high frequencies.
for subj=allsubs
    erf_osc_analysis_pow(subj, false);
end

%% Source contrast
for subj=allsubs
    erfosc_execute_pipeline('erfosc_script_jm', subj, {'dodics_gamma', true}, {'dodics_lowfreq', true});
end

%% Correlation with RT
% computing correlations and partial correlations between jitter, gamma
% power and rt.
erf_osc_analysis_corr([],[],'gamma_rt')

% correlation of rt with gamma in 3D sourcespace
for subj=allsubs
    erfosc_execute_pipeline('erfosc_script_jm', subj, {'docorr_gamma_rt', true});
end

%% Correlation Gamma ERF
for subj=allsubs
    erfosc_execute_pipeline('erfosc_script_jm', subj, {'docorrpow_lcmv',1}, {'dofreq_short', 1});
    erfosc_execute_pipeline('erfosc_script_jm', subj, {'docorrpow_lcmv_lowfreq',1}, {'dofreq_short_lowfreq', 1});
end

% do group statistics
erfosc_execute_pipeline('erfosc_script_jm', 1, {'dostat_pow_erf',1}, {'GA', 1}, {'whichFreq', 1});
erfosc_execute_pipeline('erfosc_script_jm', 1, {'dostat_pow_erf',1}, {'GA', 1}, {'whichFreq', 2});


%% Correlation ERF RT
erfosc_execute_pipeline('erfosc_script_jm', 1, {dostat_erf_rt, 'true'});




%% Plotting
plot_figures_chapter






















