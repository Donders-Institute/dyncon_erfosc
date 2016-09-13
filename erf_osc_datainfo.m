%% datainfo ERF Oscillation
% subject specific information for analysis
subjects=[];

%% pilot

pilotsubjects(1).channels = [];
pilotsubjects(1).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/01/pilot1_1200hz_20160812_01.ds';
pilotsubjects(1).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_01.mat';
pilotsubjects(1).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_01_polhemus.mat';
pilotsubjects(1).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/01/98944/ELESTJ_20141107_S32.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0001.0001.2014.11.07.16.51.27.247752.49775718.IMA';
pilotsubjects(1).icacomp = [1, 2, 8, 28];

pilotsubjects(2).channels = [];
pilotsubjects(2).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/02/pilot2_1200hz_20160815_01.ds';
pilotsubjects(2).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_02.mat';
pilotsubjects(2).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_02_polhemus.mat';
pilotsubjects(2).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/02/105520/LIEVLIE_20150616_S92.MR.LIEVLIE_FEATURE_EXP.0002.0001.2015.06.16.16.13.03.531642.240004270.IMA';
pilotsubjects(2).icacomp = [4, 7, 33];


%% Experiment

subjects(1).channels = [];
subjects(1).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/experiment/01/subj01_1200hz_20160812_01.ds';
subjects(1).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/experiment/subj_01.mat';
subjects(1).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/experiment/subj_01_polhemus.mat';
subjects(1).mri = '';
subjects(1).icacomp = [];
