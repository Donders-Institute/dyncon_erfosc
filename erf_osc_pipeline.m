erf_osc_datainfo;

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

erf_osc_analysis_tfa(subj, isPilot, 'onset');
erf_osc_analysis_tfa(subj, isPilot, 'reversal');

%% Gamma power
% estimate gamma peak frequency; create virtual gamma channel; calculate
% gamma power and gamma phase pre stimulus reversal.

erf_osc_analysis_gamma_pow(subj, isPilot);
erf_osc_analysis_gamma_freq(subj, isPilot);
erf_osc_analysis_gamma_virtualchan(subj, isPilot);



%% GLM erf

erf_osc_analysis_glm_gamma_time(subj, isPilot); % gamma power estimated at virtual gamma channel as regressor over timelocked data.

% linear regression with p1 amplitude as regressor over time-frequency
% data.
erf_osc_analysis_glm_erf_tfch(subj, isPilot, 'onset', 'onset'); % timelocked to stimulus onset, p1 amplitude of stimulus onset
erf_osc_analysis_glm_erf_tfch(subj, isPilot, 'onset', 'reversal'); % timelocked to stimulus onset, p1 amplitude of stimulus reversal
erf_osc_analysis_glm_erf_tfch(subj, isPilot, 'reversal', 'onset'); % timelocked to stimulus reversal, p1 amplitude of stimulus onset
erf_osc_analysis_glm_erf_tfch(subj, isPilot, 'reversal', 'reversal'); % timelocked to stimulus reversal, p1 amplitude of stimulus reversal
% erf_osc_analysis_erf(subj, isPilot);


% erf_osc_analysis_erf_dss_aseo(subj, isPilot)

%% statistics
% erf_osc_analysis_corr(subj, isPilot)
erf_osc_analysis_stat_glm_gamma
erf_osc_analysis_stat_tfr('high', 'onset');
erf_osc_analysis_stat_tfr('high', 'reversal');
erf_osc_analysis_stat_glm_erf('high', 'onset', 'onset');
erf_osc_analysis_stat_glm_erf('high', 'onset', 'reversal');
erf_osc_analysis_stat_glm_erf('high', 'reversal', 'onset');
erf_osc_analysis_stat_glm_erf('high', 'reversal', 'reversal');

