%% erf_osc_pipeline
% this script contains the full pipeline for the analysis of the ERF-
% oscillation experiment.
subj = 1;
isPilot = false;

%% preprocessing
% preproces MEG data: remove artifacts, apply filters etc.
% preproces MRI data: realign MRI and segment.


% interactive function: one has to inspect trial/channel variability, jump,
% muscle and blink artifacts and ICA components
existArtifact = false; % have the artifacts already been defined?
erf_osc_preprocessing_artifact(subj, isPilot, existArtifact)


% pipeline in itself, using freesurfer, workbench and realigning with
% different coordinate systems
erf_osc_preprocessing_mri

%% TFA

erf_osc_analysis_tfa(subj, isPilot);

%% Gamma power
% estimate gamma peak frequency; create virtual gamma channel; calculate
% gamma power and gamma phase pre stimulus reversal.

erf_osc_analysis_gamma_freq(subj, isPilot);

erf_osc_analysis_gamma_virtualchan(subj, isPilot)

erf_osc_analysis_gamma_pow(subj, isPilot);

% erf_osc_analysis_gamma_phase(subj, isPilot);


%% ERF
erf_osc_analysis_erf(subj, isPilot);


erf_osc_analysis_erf_dss_aseo(subj, isPilot)

%% Corelation ERF - gamma power

erf_osc_analysis_corr(subj, isPilot)
