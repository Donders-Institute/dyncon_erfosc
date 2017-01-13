%% datainfo ERF Oscillation
% subject specific information for analysis
pilotsubjects = [];
subjects=[];

%% pilot

pilotsubjects(1).channels      = [];
pilotsubjects(1).dataset       = '/project/3011085.02/raw/pilot-01/ses-meg01/pilot1_1200hz_20160812_01.ds';
pilotsubjects(1).logfile       = '/project/3011085.02/raw/pilot-01/ses-beh01/pilot_01.mat';
pilotsubjects(1).headshape     = '/project/3011085.02/raw/pilot-01/ses-beh01/pilot_01_polhemus.pos';
pilotsubjects(1).mri           = '/project/3011085.02/raw/pilot-01/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/ELESTJ_20141107_S32.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2014.11.07.16.51.27.247752.49776744.IMA';
pilotsubjects(1).segmentedmri  = '/project/3011085.02/clean/pilot-01/ses-mri01/segmentedmri';
pilotsubjects(1).ecgcomp       = [2, 8];
pilotsubjects(1).eyecomp       = [1, 15, 16, 28];
pilotsubjects(1).aseo          = [];

pilotsubjects(2).channels      = [];
pilotsubjects(2).dataset       = '/project/3011085.02/raw/pilot-02/ses-meg01/pilot2_1200hz_20160815_01.ds';
pilotsubjects(2).logfile       = '/project/3011085.02/raw/pilot-02/ses-01/pilot_02.mat';
pilotsubjects(2).headshape     = '/project/3011085.02/raw/pilot-02/ses-01/pilot_02_polhemus.pos';
pilotsubjects(2).mri           = '/project/3011085.02/raw/pilot-02/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/LIEVLIE_20150616_S92.MR.LIEVLIE_FEATURE_EXP.0002.0001.2015.06.16.16.13.03.531642.240004270.IMA';
pilotsubjects(2).segmentedmri  = '/project/3011085.02/clean/pilot-02/ses-mri01/segmentedmri';
pilotsubjects(2).ecgcomp       = [7, 16, 33];
pilotsubjects(2).eyecomp       = [4, 25];
pilotsubjects(2).aseo          = [149 ];

pilotsubjects(3).channels      = [];
pilotsubjects(3).dataset       = '/project/3011085.02/raw/pilot-03/ses-meg01/pilot3_1200hz_20160920_02.ds';
pilotsubjects(3).logfile       = '/project/3011085.02/raw/pilot-03/ses-beh01/pilot_03.mat';
pilotsubjects(3).headshape     = '/project/3011085.02/raw/pilot-03/ses-beh01/pilot_03_polhemus.pos';
pilotsubjects(3).mri           = '/project/3011085.02/raw/pilot-03/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(3).segmentedmri  = '/project/3011085.02/clean/pilot-03/ses-mri01/segmentedmri';
pilotsubjects(3).ecgcomp       = [3, 13];
pilotsubjects(3).eyecomp       = [1];
pilotsubjects(3).aseo          = [128 157;158 216];

%gamma experiment Stan van Pelt
pilotsubjects(103).channels    = [];
pilotsubjects(103).dataset     = '/project/3011085.02/raw/pilot-03/ses-meg02/pilot3_1200hz_20160920_01.ds';
pilotsubjects(103).logfile     = '/project/3011085.02/raw/pilot-03/ses-beh01/pilot_03.mat';
pilotsubjects(103).headshape   = '/project/3011085.02/raw/pilot-03/ses-beh01/pilot_03_polhemus.mat';
pilotsubjects(103).mri         = '/project/3011085.02/raw/pilot-03/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(103).segmentedmri= '/project/3011085.02/clean/pilot-03/ses-mri01/segmentedmri';
pilotsubjects(103).ecgcomp     = [3, 26];
pilotsubjects(103).eyecomp     = [1];
pilotsubjects(103).aseo        = [119 160;161 226;230 459];

pilotsubjects(4).channels      = [];
pilotsubjects(4).dataset       = '/project/3011085.02/raw/pilot-04/ses-meg01/pilot4_1200hz_20161014_01.ds';
pilotsubjects(4).logfile       = '/project/3011085.02/raw/pilot-04/ses-beh01/pilot_04.mat';
pilotsubjects(4).headshape     = '/project/3011085.02/raw/pilot-04/ses-beh01/pilot_04_polhemus.pos';
pilotsubjects(4).mri           = '/project/3011085.02/raw/pilot-04/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/JESASK_20161014_PILOT4.MR.DCCN_SKYRA.0002.0001.2016.10.14.11.32.58.155107.1286489168.IMA';
pilotsubjects(4).segmentedmri  = '/project/3011085.02/clean/pilot-04/ses-mri01/segmentedmri';
pilotsubjects(4).ecgcomp       = [3,4];
pilotsubjects(4).eyecomp       = [1];
pilotsubjects(4).aseo          = [119 161; 162 226; 230 459];


%% Experiment

subjects(1).channels           = [];
subjects(1).dataset            = '/project/3011085.02/raw/sub-01/ses-meg01/sub01ses01_3011085.02_20161130_01.ds';
subjects(1).logfile            = '/project/3011085.02/raw/sub-01/ses-beh01/logfile.mat';
subjects(1).headshape          = '/project/3011085.02/raw/sub-01/ses-beh01/polhemus.mat';
subjects(1).mri                = '/project/3011085.02/raw/sub-01/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20150915_S22.MR.DCCN_PRISMA.0002.0001.2015.09.15.11.23.24.596922.1311820.IMA';
subjects(1).segmentedmri       = '/project/3011085.02/clean/sub-01/ses-mri01/';
subjects(1).ecgcomp            = [];
subjects(1).eyecomp            = [];

subject(2).channels           = [];
subject(2).dataset            = '/project/3011085.02/raw/sub-02/ses-meg01/sub02ses01_3011085.02_20161130_01.ds';
subject(2).logfile            = '/project/3011085.02/raw/sub-02/ses-beh01/logfile.mat';
subject(2).headshape          = '/project/3011085.02/raw/sub-02/ses-beh01/polhemus.mat';
subject(2).mri                = '/project/3011085.02/raw/sub-02/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20161201_SUB02.MR.DCCN_PRISMAFIT.0002.0001.2016.12.01.13.26.11.582491.572385364.IMA';
subject(2).segmentedmri       = '/project/3011085.02/clean/sub-02/ses-mri01/';
subject(2).ecgcomp            = [];
subject(2).eyecomp            = [];

subject(3).channels           = [];
subject(3).dataset            = '/project/3011085.02/raw/sub-03/ses-meg01/sub03ses01_3011085.02_20170104_01.ds';
subject(3).logfile            = '/project/3011085.02/raw/sub-03/ses-beh01/logfile.mat';
subject(3).headshape          = '/project/3011085.02/raw/sub-03/ses-beh01/polhemus.mat';
subject(3).mri                = '/project/3011085.02/raw/sub-03/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(3).segmentedmri       = '/project/3011085.02/clean/sub-03/ses-mri01/';
subject(3).ecgcomp            = [];
subject(3).eyecomp            = [];

subject(4).channels           = [];
subject(4).dataset            = '/project/3011085.02/raw/sub-04/ses-meg01/sub04ses01_3011085.02_20170104_01.ds';
subject(4).logfile            = '/project/3011085.02/raw/sub-04/ses-beh01/logfile.mat';
subject(4).headshape          = '/project/3011085.02/raw/sub-04/ses-beh01/polhemus.mat';
subject(4).mri                = '/project/3011085.02/raw/sub-04/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/P027_t1.nii.gz';
subject(4).segmentedmri       = '/project/3011085.02/clean/sub-04/ses-mri01/';
subject(4).ecgcomp            = [];
subject(4).eyecomp            = [];

subject(5).channels           = [];
subject(5).dataset            = '/project/3011085.02/raw/sub-05/ses-meg01/sub05ses01_3011085.02_20170104_01.ds';
subject(5).logfile            = '/project/3011085.02/raw/sub-05/ses-beh01/logfile.mat';
subject(5).headshape          = '/project/3011085.02/raw/sub-05/ses-beh01/polhemus.mat';
subject(5).mri                = '/project/3011085.02/raw/sub-05/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20160517_SXX.MR.DCCN_SKYRA.0002.0001.2016.05.17.15.14.58.315238.1202276508.IMA';
subject(5).segmentedmri       = '/project/3011085.02/clean/sub-05/ses-mri01/';
subject(5).ecgcomp            = [];
subject(5).eyecomp            = [];

subject(6).channels           = [];
subject(6).dataset            = '/project/3011085.02/raw/sub-06/ses-meg01/sub06ses01_3011085.02_20170104_01.ds';
subject(6).logfile            = '/project/3011085.02/raw/sub-06/ses-beh01/logfile.mat';
subject(6).headshape          = '/project/3011085.02/raw/sub-06/ses-beh01/polhemus.mat';
subject(6).mri                = '/project/3011085.02/raw/sub-06/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/301202610_TOBNAV_S07.MR.TOBNAV_SKYRA.0006.0001.2016.11.24.12.29.16.835954.1313779829';
subject(6).segmentedmri       = '/project/3011085.02/clean/sub-06/ses-mri01/';
subject(6).ecgcomp            = [];
subject(6).eyecomp            = [];

subject(7).channels           = [];
subject(7).dataset            = '/project/3011085.02/raw/sub-07/ses-meg01/sub07ses01_3011085.02_20170104_01.ds';
subject(7).logfile            = '/project/3011085.02/raw/sub-07/ses-beh01/logfile.mat';
subject(7).headshape          = '/project/3011085.02/raw/sub-07/ses-beh01/polhemus.mat';
subject(7).mri                = '/project/3011085.02/raw/sub-07/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(7).segmentedmri       = '/project/3011085.02/clean/sub-07/ses-mri01/';
subject(7).ecgcomp            = [];
subject(7).eyecomp            = [];

subject(8).channels           = [];
subject(8).dataset            = '/project/3011085.02/raw/sub-08/ses-meg01/sub08ses01_3011085.02_20170104_01.ds';
subject(8).logfile            = '/project/3011085.02/raw/sub-08/ses-beh01/logfile.mat';
subject(8).headshape          = '/project/3011085.02/raw/sub-08/ses-beh01/polhemus.mat';
subject(8).mri                = '/project/3011085.02/raw/sub-08/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/301202610_TOBNAV_S05.MR.TOBNAV_SKYRA.0002.0001.2016.10.31.20.10.53.378701.1299671212';
subject(8).segmentedmri       = '/project/3011085.02/clean/sub-08/ses-mri01/';
subject(8).ecgcomp            = [];
subject(8).eyecomp            = [];

subject(9).channels           = [];
subject(9).dataset            = '/project/3011085.02/raw/sub-09/ses-meg01/sub09ses01_3011085.02_20170105_01.ds';
subject(9).logfile            = '/project/3011085.02/raw/sub-09/ses-beh01/logfile.mat';
subject(9).headshape          = '/project/3011085.02/raw/sub-09/ses-beh01/polhemus.mat';
subject(9).mri                = '/project/3011085.02/raw/sub-09/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170105_S09.MR.DCCN_SKYRA.0002.0001.2017.01.05.11.21.07.149215.1339150432';
subject(9).segmentedmri       = '/project/3011085.02/clean/sub-09/ses-mri01/';
subject(9).ecgcomp            = [];
subject(9).eyecomp            = [];

subject(10).channels           = [];
subject(10).dataset            = '/project/3011085.02/raw/sub-10/ses-meg01/sub10ses01_3011085.02_20170105_01.ds';
subject(10).logfile            = '/project/3011085.02/raw/sub-10/ses-beh01/logfile.mat';
subject(10).headshape          = '/project/3011085.02/raw/sub-10/ses-beh01/polhemus.mat';
subject(10).mri                = '/project/3011085.02/raw/sub-10/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170105_S10.MR.DCCN_SKYRA.0002.0001.2017.01.05.12.31.58.785395.1339167846';
subject(10).segmentedmri       = '/project/3011085.02/clean/sub-10/ses-mri01/';
subject(10).ecgcomp            = [];
subject(10).eyecomp            = [];

subject(11).channels           = [];
subject(11).dataset            = '/project/3011085.02/raw/sub-11/ses-meg01/sub11ses01_3011085.02_20170109_01.ds';
subject(11).logfile            = '/project/3011085.02/raw/sub-11/ses-beh01/logfile.mat';
subject(11).headshape          = '/project/3011085.02/raw/sub-11/ses-beh01/polhemus.mat';
subject(11).mri                = '/project/3011085.02/raw/sub-11/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/SAMCHO_20160616_COLOR20.MR.DCCN_PRISMA.0002.0001.2016.06.16.11.36.21.938505.222768478';
subject(11).segmentedmri       = '/project/3011085.02/clean/sub-11/ses-mri01/';
subject(11).ecgcomp            = [];
subject(11).eyecomp            = [];

subject(12).channels           = [];
subject(12).dataset            = '/project/3011085.02/raw/sub-12/ses-meg01/sub12ses01_3011085.02_20170111_01.ds';
subject(12).logfile            = '/project/3011085.02/raw/sub-12/ses-beh01/logfile.mat';
subject(12).headshape          = '/project/3011085.02/raw/sub-12/ses-beh01/polhemus.mat';
subject(12).mri                = '/project/3011085.02/raw/sub-12/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/BOBBRA_20160906_S29_AAT.MR.ROGMAR_PRISMA.0015.0001.2016.09.06.14.57.34.702500.352874660';
subject(12).segmentedmri       = '/project/3011085.02/clean/sub-12/ses-mri01/';
subject(12).ecgcomp            = [];
subject(12).eyecomp            = [];

subject(13).channels           = [];
subject(13).dataset            = '/project/3011085.02/raw/sub-13/ses-meg01/sub13ses01_3011085.02_20170111_01.ds';
subject(13).logfile            = '/project/3011085.02/raw/sub-13/ses-beh01/logfile.mat';
subject(13).headshape          = '/project/3011085.02/raw/sub-13/ses-beh01/polhemus.mat';
subject(13).mri                = '/project/3011085.02/raw/sub-13/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170111_S13.MR.DCCN_SKYRA.0002.0001.2017.01.11.11.51.10.945841.1344162904';
subject(13).segmentedmri       = '/project/3011085.02/clean/sub-13/ses-mri01/';
subject(13).ecgcomp            = [];
subject(13).eyecomp            = [];

subject(14).channels           = [];
subject(14).dataset            = '/project/3011085.02/raw/sub-14/ses-meg01/sub14ses01_3011085.02_20170111_01.ds';
subject(14).logfile            = '/project/3011085.02/raw/sub-14/ses-beh01/logfile.mat';
subject(14).headshape          = '/project/3011085.02/raw/sub-14/ses-beh01/polhemus.mat';
subject(14).mri                = '/project/3011085.02/raw/sub-14/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011112541298968943021';
subject(14).segmentedmri       = '/project/3011085.02/clean/sub-14/ses-mri01/';
subject(14).ecgcomp            = [];
subject(14).eyecomp            = [];

subject(15).channels           = [];
subject(15).dataset            = '/project/3011085.02/raw/sub-15/ses-meg01/sub15ses01_3011085.02_20170111_01.ds';
subject(15).logfile            = '/project/3011085.02/raw/sub-15/ses-beh01/logfile.mat';
subject(15).headshape          = '/project/3011085.02/raw/sub-15/ses-beh01/polhemus.mat';
subject(15).mri                = '/project/3011085.02/raw/sub-15/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011113580594143847289';
subject(15).segmentedmri       = '/project/3011085.02/clean/sub-15/ses-mri01/';
subject(15).ecgcomp            = [];
subject(15).eyecomp            = [];

subject(16).channels           = [];
subject(16).dataset            = '/project/3011085.02/raw/sub-16/ses-meg01/sub16ses01_3011085.02_20170112_01.ds';
subject(16).logfile            = '/project/3011085.02/raw/sub-16/ses-beh01/logfile.mat';
subject(16).headshape          = '/project/3011085.02/raw/sub-16/ses-beh01/polhemus.mat';
subject(16).mri                = '/project/3011085.02/raw/sub-16/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011211272151943534752';
subject(16).segmentedmri       = '/project/3011085.02/clean/sub-16/ses-mri01/';
subject(16).ecgcomp            = [];
subject(16).eyecomp            = [];

subject(17).channels           = [];
subject(17).dataset            = '/project/3011085.02/raw/sub-17/ses-meg01/sub17ses01_3011085.02_2017_01.ds';
subject(17).logfile            = '/project/3011085.02/raw/sub-17/ses-beh01/logfile.mat';
subject(17).headshape          = '/project/3011085.02/raw/sub-17/ses-beh01/polhemus.mat';
subject(17).mri                = '/project/3011085.02/raw/sub-17/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(17).segmentedmri       = '/project/3011085.02/clean/sub-17/ses-mri01/';
subject(17).ecgcomp            = [];
subject(17).eyecomp            = [];

subject(18).channels           = [];
subject(18).dataset            = '/project/3011085.02/raw/sub-18/ses-meg01/sub18ses01_3011085.02_20170112_01.ds';
subject(18).logfile            = '/project/3011085.02/raw/sub-18/ses-beh01/logfile.mat';
subject(18).headshape          = '/project/3011085.02/raw/sub-18/ses-beh01/polhemus.mat';
subject(18).mri                = '/project/3011085.02/raw/sub-18/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.201701121317291121285204';
subject(18).segmentedmri       = '/project/3011085.02/clean/sub-18/ses-mri01/';
subject(18).ecgcomp            = [];
subject(18).eyecomp            = [];

subject(19).channels           = [];
subject(19).dataset            = '/project/3011085.02/raw/sub-19/ses-meg01/sub19ses01_3011085.02_20170116_01.ds';
subject(19).logfile            = '/project/3011085.02/raw/sub-19/ses-beh01/logfile.mat';
subject(19).headshape          = '/project/3011085.02/raw/sub-19/ses-beh01/polhemus.mat';
subject(19).mri                = '/project/3011085.02/raw/sub-19/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(19).segmentedmri       = '/project/3011085.02/clean/sub-19/ses-mri01/';
subject(19).ecgcomp            = [];
subject(19).eyecomp            = [];

subject(20).channels           = [];
subject(20).dataset            = '/project/3011085.02/raw/sub-20/ses-meg01/sub20ses01_3011085.02_20170117_01.ds';
subject(20).logfile            = '/project/3011085.02/raw/sub-20/ses-beh01/logfile.mat';
subject(20).headshape          = '/project/3011085.02/raw/sub-20/ses-beh01/polhemus.mat';
subject(20).mri                = '/project/3011085.02/raw/sub-20/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(20).segmentedmri       = '/project/3011085.02/clean/sub-20/ses-mri01/';
subject(20).ecgcomp            = [];
subject(20).eyecomp            = [];

subject(21).channels           = [];
subject(21).dataset            = '/project/3011085.02/raw/sub-21/ses-meg01/sub21ses01_3011085.02_20170117_01.ds';
subject(21).logfile            = '/project/3011085.02/raw/sub-21/ses-beh01/logfile.mat';
subject(21).headshape          = '/project/3011085.02/raw/sub-21/ses-beh01/polhemus.mat';
subject(21).mri                = '/project/3011085.02/raw/sub-21/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(21).segmentedmri       = '/project/3011085.02/clean/sub-21/ses-mri01/';
subject(21).ecgcomp            = [];
subject(21).eyecomp            = [];

subject(22).channels           = [];
subject(22).dataset            = '/project/3011085.02/raw/sub-22/ses-meg01/sub22ses01_3011085.02_20170118_01.ds';
subject(22).logfile            = '/project/3011085.02/raw/sub-22/ses-beh01/logfile.mat';
subject(22).headshape          = '/project/3011085.02/raw/sub-22/ses-beh01/polhemus.mat';
subject(22).mri                = '/project/3011085.02/raw/sub-22/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(22).segmentedmri       = '/project/3011085.02/clean/sub-22/ses-mri01/';
subject(22).ecgcomp            = [];
subject(22).eyecomp            = [];

subject(23).channels           = [];
subject(23).dataset            = '/project/3011085.02/raw/sub-23/ses-meg01/sub23ses01_3011085.02_20170118_01.ds';
subject(23).logfile            = '/project/3011085.02/raw/sub-23/ses-beh01/logfile.mat';
subject(23).headshape          = '/project/3011085.02/raw/sub-23/ses-beh01/polhemus.mat';
subject(23).mri                = '/project/3011085.02/raw/sub-23/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(23).segmentedmri       = '/project/3011085.02/clean/sub-23/ses-mri01/';
subject(23).ecgcomp            = [];
subject(23).eyecomp            = [];

subject(24).channels           = [];
subject(24).dataset            = '/project/3011085.02/raw/sub-24/ses-meg01/sub24ses01_3011085.02_20170118_01.ds';
subject(24).logfile            = '/project/3011085.02/raw/sub-24/ses-beh01/logfile.mat';
subject(24).headshape          = '/project/3011085.02/raw/sub-24/ses-beh01/polhemus.mat';
subject(24).mri                = '/project/3011085.02/raw/sub-24/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(24).segmentedmri       = '/project/3011085.02/clean/sub-24/ses-mri01/';
subject(24).ecgcomp            = [];
subject(24).eyecomp            = [];

subject(25).channels           = [];
subject(25).dataset            = '/project/3011085.02/raw/sub-25/ses-meg01/sub25ses01_3011085.02_20170118_01.ds';
subject(25).logfile            = '/project/3011085.02/raw/sub-25/ses-beh01/logfile.mat';
subject(25).headshape          = '/project/3011085.02/raw/sub-25/ses-beh01/polhemus.mat';
subject(25).mri                = '/project/3011085.02/raw/sub-25/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(25).segmentedmri       = '/project/3011085.02/clean/sub-25/ses-mri01/';
subject(25).ecgcomp            = [];
subject(25).eyecomp            = [];

subject(26).channels           = [];
subject(26).dataset            = '/project/3011085.02/raw/sub-26/ses-meg01/sub26ses01_3011085.02_20170118_01.ds';
subject(26).logfile            = '/project/3011085.02/raw/sub-26/ses-beh01/logfile.mat';
subject(26).headshape          = '/project/3011085.02/raw/sub-26/ses-beh01/polhemus.mat';
subject(26).mri                = '/project/3011085.02/raw/sub-26/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(26).segmentedmri       = '/project/3011085.02/clean/sub-26/ses-mri01/';
subject(26).ecgcomp            = [];
subject(26).eyecomp            = [];

subject(27).channels           = [];
subject(27).dataset            = '/project/3011085.02/raw/sub-27/ses-meg01/sub27ses01_3011085.02_20170123_01.ds';
subject(27).logfile            = '/project/3011085.02/raw/sub-27/ses-beh01/logfile.mat';
subject(27).headshape          = '/project/3011085.02/raw/sub-27/ses-beh01/polhemus.mat';
subject(27).mri                = '/project/3011085.02/raw/sub-27/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(27).segmentedmri       = '/project/3011085.02/clean/sub-27/ses-mri01/';
subject(27).ecgcomp            = [];
subject(27).eyecomp            = [];

subject(28).channels           = [];
subject(28).dataset            = '/project/3011085.02/raw/sub-28/ses-meg01/sub28ses01_3011085.02_20170125_01.ds';
subject(28).logfile            = '/project/3011085.02/raw/sub-28/ses-beh01/logfile.mat';
subject(28).headshape          = '/project/3011085.02/raw/sub-28/ses-beh01/polhemus.mat';
subject(28).mri                = '/project/3011085.02/raw/sub-28/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(28).segmentedmri       = '/project/3011085.02/clean/sub-28/ses-mri01/';
subject(28).ecgcomp            = [];
subject(28).eyecomp            = [];

subject(29).channels           = [];
subject(29).dataset            = '/project/3011085.02/raw/sub-29/ses-meg01/sub29ses01_3011085.02_20170125_01.ds';
subject(29).logfile            = '/project/3011085.02/raw/sub-29/ses-beh01/logfile.mat';
subject(29).headshape          = '/project/3011085.02/raw/sub-29/ses-beh01/polhemus.mat';
subject(29).mri                = '/project/3011085.02/raw/sub-29/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(29).segmentedmri       = '/project/3011085.02/clean/sub-29/ses-mri01/';
subject(29).ecgcomp            = [];
subject(29).eyecomp            = [];

subject(30).channels           = [];
subject(30).dataset            = '/project/3011085.02/raw/sub-30/ses-meg01/sub30ses01_3011085.02_20170125_01.ds';
subject(30).logfile            = '/project/3011085.02/raw/sub-30/ses-beh01/logfile.mat';
subject(30).headshape          = '/project/3011085.02/raw/sub-30/ses-beh01/polhemus.mat';
subject(30).mri                = '/project/3011085.02/raw/sub-30/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(30).segmentedmri       = '/project/3011085.02/clean/sub-30/ses-mri01/';
subject(30).ecgcomp            = [];
subject(30).eyecomp            = [];

subject(31).channels           = [];
subject(31).dataset            = '/project/3011085.02/raw/sub-31/ses-meg01/sub31ses01_3011085.02_20170126_01.ds';
subject(31).logfile            = '/project/3011085.02/raw/sub-31/ses-beh01/logfile.mat';
subject(31).headshape          = '/project/3011085.02/raw/sub-31/ses-beh01/polhemus.mat';
subject(31).mri                = '/project/3011085.02/raw/sub-31/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(31).segmentedmri       = '/project/3011085.02/clean/sub-31/ses-mri01/';
subject(31).ecgcomp            = [];
subject(31).eyecomp            = [];

subject(32).channels           = [];
subject(32).dataset            = '/project/3011085.02/raw/sub-32/ses-meg01/sub32ses01_3011085.02_20170126_01.ds';
subject(32).logfile            = '/project/3011085.02/raw/sub-32/ses-beh01/logfile.mat';
subject(32).headshape          = '/project/3011085.02/raw/sub-32/ses-beh01/polhemus.mat';
subject(32).mri                = '/project/3011085.02/raw/sub-32/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(32).segmentedmri       = '/project/3011085.02/clean/sub-32/ses-mri01/';
subject(32).ecgcomp            = [];
subject(32).eyecomp            = [];

subject(33).channels           = [];
subject(33).dataset            = '/project/3011085.02/raw/sub-33/ses-meg01/sub33ses01_3011085.02_20170126_01.ds';
subject(33).logfile            = '/project/3011085.02/raw/sub-33/ses-beh01/logfile.mat';
subject(33).headshape          = '/project/3011085.02/raw/sub-33/ses-beh01/polhemus.mat';
subject(33).mri                = '/project/3011085.02/raw/sub-33/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/';
subject(33).segmentedmri       = '/project/3011085.02/clean/sub-33/ses-mri01/';
subject(33).ecgcomp            = [];
subject(33).eyecomp            = [];