%% datainfo ERF Oscillation
% subject specific information for analysis
pilotsubjects = [];
subjects=[];

%% pilot

pilotsubjects(1).channels = [];
pilotsubjects(1).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/01/pilot1_1200hz_20160812_01.ds';
pilotsubjects(1).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_01.mat';
pilotsubjects(1).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_01_polhemus.pos';
pilotsubjects(1).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/01/98944/ELESTJ_20141107_S32.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0001.0001.2014.11.07.16.51.27.247752.49775718.IMA';
pilotsubjects(1).segmentedmri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/01/segmentedmri.mat';
pilotsubjects(1).ecgcomp = [2, 8];
pilotsubjects(1).eyecomp = [1, 15, 16, 28];

pilotsubjects(2).channels = [];
pilotsubjects(2).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/02/pilot2_1200hz_20160815_01.ds';
pilotsubjects(2).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_02.mat';
pilotsubjects(2).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_02_polhemus.pos';
pilotsubjects(2).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/02/105520/LIEVLIE_20150616_S92.MR.LIEVLIE_FEATURE_EXP.0002.0001.2015.06.16.16.13.03.531642.240004270.IMA';
pilotsubjects(2).segmentedmri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/02/segmentedmri.mat';
pilotsubjects(2).ecgcomp = [7, 16, 33];
pilotsubjects(2).eyecomp = [4, 25];

pilotsubjects(3).channels = [];
pilotsubjects(3).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/03/pilot3_1200hz_20160920_02.ds';
pilotsubjects(3).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_03.mat';
pilotsubjects(3).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_03_polhemus.pos';
pilotsubjects(3).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/03/101497/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(3).segmentedmri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/03/segmentedmri';
pilotsubjects(3).ecgcomp = [3, 13];
pilotsubjects(3).eyecomp = [1];

%gamma experiment Stan van Pelt
pilotsubjects(103).channels = [];
pilotsubjects(103).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/03/pilot3_1200hz_20160920_01.ds';
pilotsubjects(103).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_03.mat';
pilotsubjects(103).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_03_polhemus.mat';
pilotsubjects(103).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/03/101497/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(103).segmentedmri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/03/segmentedmri';
pilotsubjects(103).ecgcomp = [3, 26];
pilotsubjects(103).eyecomp = [1];

pilotsubjects(4).channels = [];
pilotsubjects(4).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot/04/pilot4_1200hz_20161014_01.ds';
pilotsubjects(4).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_04.mat';
pilotsubjects(4).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot/pilot_04_polhemus.pos';
pilotsubjects(4).mri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/04/20161014_pilot4/JESASK_20161014_PILOT4.MR.DCCN_SKYRA.0002.0001.2016.10.14.11.32.58.155107.1286489168.IMA';
pilotsubjects(4).segmentedmri = '/home/electromag/matves/Data/ERF_oscillation/sMRI/pilot/04/segmentedmri';
pilotsubjects(4).ecgcomp = [3,4];
pilotsubjects(4).eyecomp = [1];


%% Experiment

subjects(1).channels = [];
subjects(1).dataset = '/home/electromag/matves/Data/ERF_oscillation/raw_data/experiment/01/subj01_1200hz_20160812_01.ds';
subjects(1).logfile = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/experiment/subj_01.mat';
subjects(1).headshape = '/home/electromag/matves/Data/ERF_oscillation/behavioral_log/experiment/subj_01_polhemus.mat';
subjects(1).mri = '';
subjects(1).ecgcomp = [];
subjects(1).eyecomp = [];

