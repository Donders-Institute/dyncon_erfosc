%% datainfo ERF Oscillation
% subject specific information for analysis
pilotsubjects = [];
subjects=[];

%% pilot

pilotsubjects(1).channels         = [];
pilotsubjects(1).dataset       = '/project/3011085.02/raw/pilot-001/ses-meg01/pilot1_1200hz_20160812_01.ds';
pilotsubjects(1).logfile       = '/project/3011085.02/raw/pilot-001/ses-beh01/pilot_01.mat';
pilotsubjects(1).headshape     = '/project/3011085.02/raw/pilot-001/polhemus_erfosc_pilot-001.pos';
pilotsubjects(1).mri           = '/project/3011085.02/raw/pilot-001/ses-mri01/002-t1_mpr_ns_PH8/ELESTJ_20141107_S32.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2014.11.07.16.51.27.247752.49776744.IMA';
pilotsubjects(1).mridir        = '/project/3011085.02/processed/pilot-001/ses-mri01/';
pilotsubjects(1).ecgcomp       = [2, 12];
pilotsubjects(1).eyecomp       = [1, 9];
pilotsubjects(1).aseo          = [];

pilotsubjects(2).channels         = [];
pilotsubjects(2).dataset       = '/project/3011085.02/raw/pilot-002/ses-meg01/pilot2_1200hz_20160815_01.ds';
pilotsubjects(2).logfile       = '/project/3011085.02/raw/pilot-002/ses-beh01/pilot_02.mat';
pilotsubjects(2).headshape     = '/project/3011085.02/raw/pilot-002/polhemus_erfosc_pilot-002.pos';
pilotsubjects(2).mri           = '/project/3011085.02/raw/pilot-002/ses-mri01/002-t1_mprage_sag_p2_iso_1.0/LIEVLIE_20150616_S92.MR.LIEVLIE_FEATURE_EXP.0002.0001.2015.06.16.16.13.03.531642.240004270.IMA';
pilotsubjects(2).mridir        = '/project/3011085.02/processed/pilot-002/ses-mri01/';
pilotsubjects(2).ecgcomp       = [7, 16, 33];
pilotsubjects(2).eyecomp       = [4, 25];
pilotsubjects(2).aseo          = [149];

pilotsubjects(3).channels         = [];
pilotsubjects(3).dataset       = '/project/3011085.02/raw/pilot-003/ses-meg01/pilot3_1200hz_20160920_02.ds';
pilotsubjects(3).logfile       = '/project/3011085.02/raw/pilot-003/ses-beh01/pilot_03.mat';
pilotsubjects(3).headshape     = '/project/3011085.02/raw/pilot-003/polhemus_erfosc_pilot-003.pos';
pilotsubjects(3).mri           = '/project/3011085.02/raw/pilot-003/ses-mri01/011-t1_mprage_sag_p2_iso_1.0/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(3).mridir        = '/project/3011085.02/processed/pilot-003/ses-mri01/';
pilotsubjects(3).ecgcomp       = [3, 13];
pilotsubjects(3).eyecomp       = [1];
pilotsubjects(3).aseo          = [128 157;158 216];

%gamma experiment Stan van Pelt
pilotsubjects(103).channels       = [];
pilotsubjects(103).dataset     = '/project/3011085.02/raw/pilot-003/ses-meg02/pilot3_1200hz_20160920_01.ds';
pilotsubjects(103).logfile     = '/project/3011085.02/raw/pilot-003/ses-beh01/pilot_03.mat';
pilotsubjects(103).headshape   = '/project/3011085.02/raw/pilot-003/polhemus_erfosc_pilot-003.mat';
pilotsubjects(103).mri         = '/project/3011085.02/raw/pilot-003/ses-mri01/011-t1_mprage_sag_p2_iso_1.0/SARFAB_07042016_KD_SUB5.MR.SARFAB_SKYRA.0011.0001.2016.04.07.18.27.54.54019.1168606464.IMA';
pilotsubjects(103).segmentedmri= '/project/3011085.02/processed/pilot-003/ses-mri01/';
pilotsubjects(103).ecgcomp     = [3, 26];
pilotsubjects(103).eyecomp     = [1];
pilotsubjects(103).aseo        = [119 160;161 226;230 459];

pilotsubjects(4).channels         = [];
pilotsubjects(4).dataset       = '/project/3011085.02/raw/pilot-004/ses-meg01/pilot4_1200hz_20161014_01.ds';
pilotsubjects(4).logfile       = '/project/3011085.02/raw/pilot-004/ses-beh01/pilot_04.mat';
pilotsubjects(4).headshape     = '/project/3011085.02/raw/pilot-004/polhemus_erfosc_pilot-004.pos';
pilotsubjects(4).mri           = '/project/3011085.02/raw/pilot-004/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/JESASK_20161014_PILOT4.MR.DCCN_SKYRA.0002.0001.2016.10.14.11.32.58.155107.1286489168.IMA';
pilotsubjects(4).mridir        = '/project/3011085.02/processed/pilot-004/ses-mri01/';
pilotsubjects(4).ecgcomp       = [3,4];
pilotsubjects(4).eyecomp       = [1];
% pilotsubjects(4).aseo          = [119 161; 162 226; 230 459];
pilotsubjects(4).aseo          = [0.048 0.083; 0.083 0.1375; 0.1375 0.4];


%% Experiment
badsubjects = [10]; %5 (dental fillings left temporal),  8 eyesaccades, 10 (corrupt logfile, only 9 blocks completed, many eye blinks)
allsubs = 1:33;
allsubs(badsubjects) = [];


subjects(1).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(1).dataset            = '/project/3011085.02/raw/sub-001/ses-meg01/sub01ses01_3011085.02_20161130_01.ds';
subjects(1).logfile            = '/project/3011085.02/raw/sub-001/ses-beh01/sub01ses01.mat';
subjects(1).headshape          = '/project/3011085.02/raw/sub-001/polhemus_erfosc_sub-001.pos';
subjects(1).mri                = '/project/3011085.02/raw/sub-001/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20150915_S22.MR.DCCN_PRISMA.0002.0001.2015.09.15.11.23.24.596922.1311820.IMA';
subjects(1).mridir             = '/project/3011085.02/processed/sub-001/ses-mri01/';
subjects(1).ecgcomp            = [5, 11];
subjects(1).eyecomp            = [3];
subjects(1).badtrials          = []; 
subjects(1).aseo               = [0.05167 0.1033; 0.1042 0.1683; 0.1858  0.5475]; 
subjects(1).gammaband          = [50 70];%TFR [60 70];powspec
subjects(1).gammapeak           = 66;

subjects(2).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(2).dataset            = '/project/3011085.02/raw/sub-002/ses-meg01/sub02ses01_3011085.02_20161130_01.ds';
subjects(2).logfile            = '/project/3011085.02/raw/sub-002/ses-beh01/sub02ses01.mat';
subjects(2).headshape          = '/project/3011085.02/raw/sub-002/polhemus_erfosc_sub-002.pos';
subjects(2).mri                = '/project/3011085.02/raw/sub-002/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20161201_SUB02.MR.DCCN_PRISMAFIT.0002.0001.2016.12.01.13.26.11.582491.572385364.IMA';
subjects(2).mridir             = '/project/3011085.02/processed/sub-002/ses-mri01/';
subjects(2).ecgcomp            = [10, 18];
subjects(2).eyecomp            = [2, 4];
subjects(2).badtrials          = [140]; % bad trials were identified with ICA and imagesc. they are removed by adding their sampleinfo data to artfctdef.eyeblink.artifact;
subjects(2).gammaband          = [56 70];%TFR [34 46; 52 70];powspec
subjects(2).gammapeak          = [40; 62]; % 62 highest

subjects(3).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(3).dataset            = '/project/3011085.02/raw/sub-003/ses-meg01/sub003ses01_3011085.02_20170220_01.ds';
subjects(3).logfile            = '/project/3011085.02/raw/sub-003/ses-beh01/sub03ses01.mat';
subjects(3).headshape          = '/project/3011085.02/raw/sub-003/polhemus_erfosc_sub-003.pos';
subjects(3).mri                = '/project/3011085.02/raw/sub-003/ses-mri01/010-t1_mprage_sag_p2_iso_1.0/STETHE_26012017_S06.MR.STETHE_PRISMAFIT.0010.0001.2017.01.26.13.29.47.463711.620607576.IMA';
subjects(3).mridir             = '/project/3011085.02/processed/sub-003/ses-mri01/';
subjects(3).ecgcomp            = [5, 7];
subjects(3).eyecomp            = [3];
subjects(3).badtrials          = [];
subjects(3).gammaband          = [48 62];% TFR [26 34; 42 68];powspec
subjects(3).gammapeak          = [30; 56]; % 56 highest

subjects(4).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(4).dataset            = '/project/3011085.02/raw/sub-004/ses-meg01/sub04ses01_3011085.02_20170104_01.ds';
subjects(4).logfile            = '/project/3011085.02/raw/sub-004/ses-beh01/sub04ses01.mat';
subjects(4).headshape          = '/project/3011085.02/raw/sub-004/polhemus_erfosc_sub-004.pos';
subjects(4).mri                = '/project/3011085.02/raw/sub-004/ses-mri01/P027_t1.nii.gz';
subjects(4).mridir             = '/project/3011085.02/processed/sub-004/ses-mri01/';
subjects(4).ecgcomp            = [6, 9];
subjects(4).eyecomp            = [2];
subjects(4).badtrials          = [];
subjects(4).gammaband          = [46 62];%TFR [40 64];powspec
subjects(4).gammapeak          = 54;

subjects(5).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(5).dataset            = '/project/3011085.02/raw/sub-005/ses-meg01/sub05ses01_3011085.02_20170104_01.ds';
subjects(5).logfile            = '/project/3011085.02/raw/sub-005/ses-beh01/sub05ses01.mat';
subjects(5).headshape          = '/project/3011085.02/raw/sub-005/polhemus_erfosc_sub-005.pos';
subjects(5).mri                = '/project/3011085.02/raw/sub-005/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20160517_SXX.MR.DCCN_SKYRA.0002.0001.2016.05.17.15.14.58.315238.1202276508.IMA';
subjects(5).mridir             = '/project/3011085.02/processed/sub-005/ses-mri01/';
subjects(5).ecgcomp            = [9];
subjects(5).eyecomp            = [6, 8]; % tooth filling left temporal, high noise levels in MLT sensors. mark components 1. CHECK FOR OTHER COMPONENTS. MANY SUSPECTED
subjects(5).badtrials          = [314, 321]; % bad trials were identified with ICA and imagesc. they are removed by adding their sampleinfo data to artfctdef.eyeblink.artifact;
subjects(5).gammaband          = [46 60];%TFR [42 60];powspec
subjects(5).gammapeak          = 54;

subjects(6).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(6).dataset            = '/project/3011085.02/raw/sub-006/ses-meg01/sub06ses01_3011085.02_20170104_01.ds';
subjects(6).logfile            = '/project/3011085.02/raw/sub-006/ses-beh01/sub06ses01.mat';
subjects(6).headshape          = '/project/3011085.02/raw/sub-006/polhemus_erfosc_sub-006.pos';
subjects(6).mri                = '/project/3011085.02/raw/sub-006/ses-mri01/007-t1_mprage_sag_p2_iso_1.0/301202610_TOBNAV_S07.MR.TOBNAV_SKYRA.0007.0001.2016.11.24.12.29.16.835954.1313789047.IMA';
subjects(6).mridir             = '/project/3011085.02/processed/sub-006/ses-mri01/';
subjects(6).ecgcomp            = [9, 34];
subjects(6).eyecomp            = [1];
subjects(6).badtrials          = [];
subjects(6).gammaband          = [42 54]; %TFR [46 62]; powspec % absent
subjects(6).gammapeak          = [46];

subjects(7).channels           = {'MEG', '-MRT32','-MLF23','HLC*'};
subjects(7).dataset            = '/project/3011085.02/raw/sub-007/ses-meg01/sub07ses01_3011085.02_20170104_01.ds';
subjects(7).logfile            = '/project/3011085.02/raw/sub-007/ses-beh01/sub07ses01.mat';
subjects(7).headshape          = '/project/3011085.02/raw/sub-007/polhemus_erfosc_sub-007.pos';
subjects(7).mri                = '/project/3011085.02/raw/sub-007/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011710180816482222404.IMA';
subjects(7).mridir             = '/project/3011085.02/processed/sub-007/ses-mri01/';
subjects(7).ecgcomp            = [9, 34];
subjects(7).eyecomp            = [1, 2, 6, 7];
subjects(7).badtrials          = [];
subjects(7).gammaband          = [50 64];%TFR [46 68];powspec
subjects(7).gammapeak           = 56;

% huge amount of eye blinks
subjects(8).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(8).dataset            = '/project/3011085.02/raw/sub-008/ses-meg01/sub08ses01_3011085.02_20170104_01.ds';
subjects(8).logfile            = '/project/3011085.02/raw/sub-008/ses-beh01/sub08ses01.mat';
subjects(8).headshape          = '/project/3011085.02/raw/sub-008/polhemus_erfosc_sub-008.pos';
subjects(8).mri                = '/project/3011085.02/raw/sub-008/ses-mri01/002-t1_mprage_sag_p2_iso_1.0/301202610_TOBNAV_S05.MR.TOBNAV_SKYRA.0002.0001.2016.10.31.20.10.53.378701.1299671212.IMA';
subjects(8).mridir             = '/project/3011085.02/processed/sub-008/ses-mri01/';
subjects(8).ecgcomp            = [8, 9];
subjects(8).eyecomp            = [1];
subjects(8).badtrials          = [];
subjects(8).gammaband          = [46 64];%TFR [46 62];powspec
subjects(8).gammapeak          = 54;

subjects(9).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(9).dataset            = '/project/3011085.02/raw/sub-009/ses-meg01/sub09ses01_3011085.02_20170105_01.ds';
subjects(9).logfile            = '/project/3011085.02/raw/sub-009/ses-beh01/sub09ses01.mat';
subjects(9).headshape          = '/project/3011085.02/raw/sub-009/polhemus_erfosc_sub-009.pos';
subjects(9).mri                = '/project/3011085.02/raw/sub-009/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170105_S09.MR.DCCN_SKYRA.0002.0001.2017.01.05.11.21.07.149215.1339150432.IMA';
subjects(9).mridir             = '/project/3011085.02/processed/sub-009/ses-mri01/';
subjects(9).ecgcomp            = [9, 34];
subjects(9).eyecomp            = [1, 2];
subjects(9).badtrials          = [];
subjects(9).gammaband          = [50 64];%TFR %[50 68];powspec
subjects(9).gammapeak          = 54;

% subject 10 bad subject: corrupt logfile, could only do 9 blocks; many
% eyeblinks. Discard subject
subjects(10).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(10).dataset            = '/project/3011085.02/raw/sub-010/ses-meg01/sub10ses01_3011085.02_20170105_01.ds';
subjects(10).logfile            = '/project/3011085.02/raw/sub-010/ses-beh01/sub10ses01.mat';
subjects(10).headshape          = '/project/3011085.02/raw/sub-010/polhemus_erfosc_sub-010.pos';
subjects(10).mri                = '/project/3011085.02/raw/sub-010/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170105_S10.MR.DCCN_SKYRA.0002.0001.2017.01.05.12.31.58.785395.1339167846.IMA';
subjects(10).mridir             = '/project/3011085.02/processed/sub-010/ses-mri01/';
subjects(10).ecgcomp            = [];
subjects(10).eyecomp            = [];
subjects(10).badtrials          = [];

subjects(11).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(11).dataset            = '/project/3011085.02/raw/sub-011/ses-meg01/sub11ses01_3011085.02_20170109_01.ds';
subjects(11).logfile            = '/project/3011085.02/raw/sub-011/ses-beh01/sub11ses01.mat';
subjects(11).headshape          = '/project/3011085.02/raw/sub-011/polhemus_erfosc_sub-011.pos';
subjects(11).mri                = '/project/3011085.02/raw/sub-011/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/SAMCHO_20160616_COLOR20.MR.DCCN_PRISMA.0002.0001.2016.06.16.11.36.21.938505.222768478.IMA';
subjects(11).mridir             = '/project/3011085.02/processed/sub-011/ses-mri01/';
subjects(11).ecgcomp            = [2, 7, 18];
subjects(11).eyecomp            = [1, 15];
subjects(11).badtrials          = [];
subjects(11).gammaband          = [50 66];%TFR [46 70];powpsec
subjects(11).gammapeak          = 56;

subjects(12).channels           = {'MEG', '-MRT32', '-MRC11', '-MZF03', 'HLC*'};
subjects(12).dataset            = '/project/3011085.02/raw/sub-012/ses-meg01/sub012ses01_3011085.02_20170222_01.ds';
subjects(12).logfile            = '/project/3011085.02/raw/sub-012/ses-beh01/sub12ses01.mat';
subjects(12).headshape          = '/project/3011085.02/raw/sub-012/polhemus_erfosc_sub-012.pos';
subjects(12).mri                = '/project/3011085.02/raw/sub-012/ses-mri01/007-t1_mprage_sag_p2_iso_1.0/301202610_TOBNAV_S15.MR.TOBNAV_SKYRA.0007.0001.2017.02.10.18.24.36.801369.1390682261.IMA';
subjects(12).mridir             = '/project/3011085.02/processed/sub-012/ses-mri01/';
subjects(12).ecgcomp            = [4, 15];
subjects(12).eyecomp            = [1,2];
subjects(12).badtrials          = [];
subjects(12).gammaband          = [42 58];%TFR [40 60];powspec
subjects(12).gammapeak          = 52;

subjects(13).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(13).dataset            = '/project/3011085.02/raw/sub-013/ses-meg01/sub13ses01_3011085.02_20170111_01.ds';
subjects(13).logfile            = '/project/3011085.02/raw/sub-013/ses-beh01/sub13ses01.mat';
subjects(13).headshape          = '/project/3011085.02/raw/sub-013/polhemus_erfosc_sub-013.pos';
subjects(13).mri                = '/project/3011085.02/raw/sub-013/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170111_S13.MR.DCCN_SKYRA.0002.0001.2017.01.11.11.51.10.945841.1344162904.IMA';
subjects(13).mridir             = '/project/3011085.02/processed/sub-013/ses-mri01/';
subjects(13).ecgcomp            = [12, 20];
subjects(13).eyecomp            = [1,2,3];
subjects(13).badtrials          = []; 
subjects(13).gammaband          = [36 50];% TFR [34 48];powspec
subjects(13).gammapeak          = 42;

subjects(14).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(14).dataset            = '/project/3011085.02/raw/sub-014/ses-meg01/sub14ses01_3011085.02_20170111_01.ds';
subjects(14).logfile            = '/project/3011085.02/raw/sub-014/ses-beh01/sub14ses01.mat';
subjects(14).headshape          = '/project/3011085.02/raw/sub-014/polhemus_erfosc_sub-014.pos';
subjects(14).mri                = '/project/3011085.02/raw/sub-014/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011112541298968943021.IMA';
subjects(14).mridir             = '/project/3011085.02/processed/sub-014/ses-mri01/';
subjects(14).ecgcomp            = [5, 8, 13];
subjects(14).eyecomp            = [1];
subjects(14).badtrials          = [];
subjects(14).gammaband          = [58 70];%TFR [54 74]; powspec
subjects(14).gammapeak          = 64;% peakFreq _gamma WRONG??

subjects(15).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(15).dataset            = '/project/3011085.02/raw/sub-015/ses-meg01/sub15ses01_3011085.02_20170111_01.ds';
subjects(15).logfile            = '/project/3011085.02/raw/sub-015/ses-beh01/sub15ses01.mat';
subjects(15).headshape          = '/project/3011085.02/raw/sub-015/polhemus_erfosc_sub-015.pos';
subjects(15).mri                = '/project/3011085.02/raw/sub-015/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011113580594143847289.IMA';
subjects(15).mridir             = '/project/3011085.02/processed/sub-015/ses-mri01/';
subjects(15).ecgcomp            = [3];
subjects(15).eyecomp            = [2];
subjects(15).badtrials          = []; 
subjects(15).gammaband          = [62 82];%TFR [64 82];powspec % absent
subjects(15).gammapeak          = [80];

subjects(16).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(16).dataset            = '/project/3011085.02/raw/sub-016/ses-meg01/sub16ses01_3011085.02_20170112_01.ds';
subjects(16).logfile            = '/project/3011085.02/raw/sub-016/ses-beh01/sub16ses01.mat';
subjects(16).headshape          = '/project/3011085.02/raw/sub-016/polhemus_erfosc_sub-016.pos';
subjects(16).mri                = '/project/3011085.02/raw/sub-016/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011211272151943534752.IMA';
subjects(16).mridir             = '/project/3011085.02/processed/sub-016/ses-mri01/';
subjects(16).ecgcomp            = [4, 6, 30];
subjects(16).eyecomp            = [2];
subjects(16).badtrials          = [];
subjects(16).gammaband          = [50 66];%TFR [46 70];powspec
subjects(16).gammapeak          = 58;

subjects(17).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(17).dataset            = '/project/3011085.02/raw/sub-017/ses-meg01/sub017ses01_3011085.02_20170222_01.ds';
subjects(17).logfile            = '/project/3011085.02/raw/sub-017/ses-beh01/sub17ses01.mat';
subjects(17).headshape          = '/project/3011085.02/raw/sub-017/polhemus_erfosc_sub-017.pos';
subjects(17).mri                = '/project/3011085.02/raw/sub-017/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017022312154313992946477.IMA';
subjects(17).mridir             = '/project/3011085.02/processed/sub-017/ses-mri01/';
subjects(17).ecgcomp            = [10, 17];
subjects(17).eyecomp            = [3];
subjects(17).badtrials          = [];
subjects(17).gammaband          = [66 78];%TFR [60 84];powspec
subjects(17).gammapeak          = 70;

subjects(18).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(18).dataset            = '/project/3011085.02/raw/sub-018/ses-meg01/sub18ses01_3011085.02_20170112_01.ds';
subjects(18).logfile            = '/project/3011085.02/raw/sub-018/ses-beh01/sub18ses01.mat';
subjects(18).headshape          = '/project/3011085.02/raw/sub-018/polhemus_erfosc_sub-018.pos';
subjects(18).mri                = '/project/3011085.02/raw/sub-018/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.201701121317291121285204.IMA';
subjects(18).mridir             = '/project/3011085.02/processed/sub-018/ses-mri01/';
subjects(18).ecgcomp            = [5, 7];
subjects(18).eyecomp            = [2, 9];
subjects(18).badtrials          = [];
subjects(18).gammaband          = [58 84]%;TFR [60 70; 74 96];powspec % absent? also in TFR
subjects(18).gammapeak          = [62; 86];

subjects(19).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(19).dataset            = '/project/3011085.02/raw/sub-019/ses-meg01/sub19ses01_3011085.02_20170116_01.ds';
subjects(19).logfile            = '/project/3011085.02/raw/sub-019/ses-beh01/sub19ses01.mat';
subjects(19).headshape          = '/project/3011085.02/raw/sub-019/polhemus_erfosc_sub-019.pos';
subjects(19).mri                = '/project/3011085.02/raw/sub-019/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017011616052260254962913.IMA';
subjects(19).mridir             = '/project/3011085.02/processed/sub-019/ses-mri01/';
subjects(19).ecgcomp            = [1, 4, 13];
subjects(19).eyecomp            = [2];
subjects(19).badtrials          = [];
subjects(19).gammaband          = [50 68];%TFR [50 70];powspec
subjects(19).gammapeak          = 58;

subjects(20).channels           = {'MEG', '-MRT32', '-MLP31', 'HLC*'};
subjects(20).dataset            = '/project/3011085.02/raw/sub-020/ses-meg01/sub20ses01_3011085.02_20170117_01.ds';
subjects(20).logfile            = '/project/3011085.02/raw/sub-020/ses-beh01/sub20ses01.mat';
subjects(20).headshape          = '/project/3011085.02/raw/sub-020/polhemus_erfosc_sub-020.pos';
subjects(20).mri                = '/project/3011085.02/raw/sub-020/ses-mri01/002-t1_mpr_ns_PH8/PAUGAA_20150826_100903.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2015.08.26.10.00.01.246391.106902594.IMA';
subjects(20).mridir             = '/project/3011085.02/processed/sub-020/ses-mri01/';
subjects(20).ecgcomp            = [9, 21];
subjects(20).eyecomp            = [2];
subjects(20).badtrials          = [299, 300]; % bad trials were identified with ICA and imagesc. they are removed by adding their sampleinfo data to artfctdef.eyeblink.artifact;
subjects(20).gammaband          = [54 68];% TFR [50 72];powspec
subjects(20).gammapeak          = 60;

subjects(21).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(21).dataset            = '/project/3011085.02/raw/sub-021/ses-meg01/sub021ses01_3011085.02_20170223_01.ds';
subjects(21).logfile            = '/project/3011085.02/raw/sub-021/ses-beh01/sub21ses01.mat';
subjects(21).headshape          = '/project/3011085.02/raw/sub-021/polhemus_erfosc_sub-021.pos';
subjects(21).mri                = '/project/3011085.02/raw/sub-021/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017022315382267121048627.IMA';
subjects(21).mridir             = '/project/3011085.02/processed/sub-021/ses-mri01/';
subjects(21).ecgcomp            = [10, 11];
subjects(21).eyecomp            = [2];
subjects(21).badtrials          = [446, 447];
subjects(21).gammaband          = [58 72];%TFR  [50 76];powspec
subjects(21).gammapeak          = 64;

subjects(22).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(22).dataset            = '/project/3011085.02/raw/sub-022/ses-meg01/sub22ses01_3011085.02_20170118_01.ds';
subjects(22).logfile            = '/project/3011085.02/raw/sub-022/ses-beh01/sub22ses01.mat';
subjects(22).headshape          = '/project/3011085.02/raw/sub-022/polhemus_erfosc_sub-022.pos';
subjects(22).mri                = '/project/3011085.02/raw/sub-022/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20170118_S22.MR.DCCN_SKYRA.0002.0001.2017.01.18.10.45.11.593116.1351346244.IMA';
subjects(22).mridir             = '/project/3011085.02/processed/sub-022/ses-mri01/';
subjects(22).ecgcomp            = [1, 2];
subjects(22).eyecomp            = [3];
subjects(22).badtrials          = [];
subjects(22).gammaband          = [54 70];% TFR [46 74];powspec
subjects(22).gammapeak          = 62;

subjects(23).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(23).dataset            = '/project/3011085.02/raw/sub-023/ses-meg01/sub23ses01_3011085.02_20170118_01.ds';
subjects(23).logfile            = '/project/3011085.02/raw/sub-023/ses-beh01/sub23ses01.mat';
subjects(23).headshape          = '/project/3011085.02/raw/sub-023/polhemus_erfosc_sub-023.pos';
subjects(23).mri                = '/project/3011085.02/raw/sub-023/ses-mri01/002-t1_mpr_ns_PH8/VITPIA_20131002_S12.MR.FCDC_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2013.10.02.13.41.11.952197.56981592.IMA';
subjects(23).mridir             = '/project/3011085.02/processed/sub-023/ses-mri01/';
subjects(23).ecgcomp            = [6, 20];
subjects(23).eyecomp            = [1, 5, 8];
subjects(23).badtrials          = [];
subjects(23).gammaband          = [60 74];%TFR  [56 74];powspec
subjects(23).gammapeak          = 66;

subjects(24).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(24).dataset            = '/project/3011085.02/raw/sub-024/ses-meg01/sub024ses01_3011085.02_20170302_01.ds';
subjects(24).logfile            = '/project/3011085.02/raw/sub-024/ses-beh01/sub24ses01.mat';
subjects(24).headshape          = '/project/3011085.02/raw/sub-024/polhemus_erfosc_sub-024.pos';
subjects(24).mri                = '/project/3011085.02/raw/sub-024/ses-mri01/016-t1_mprage_sag_p2_iso_1.0/3012027_01_LONSCH_QUANEU_S29.MR.LONSCH_PRISMAFIT.0016.0001.2016.12.08.14.00.22.549596.578717804.IMA';
subjects(24).mridir             = '/project/3011085.02/processed/sub-024/ses-mri01/';
subjects(24).ecgcomp            = [4, 5];
subjects(24).eyecomp            = [1, 2, 11];
subjects(24).badtrials          = [];
subjects(24).gammaband          = [46 62];%TFR [40 66];powspec
subjects(24).gammapeak          = 52;

subjects(25).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(25).dataset            = '/project/3011085.02/raw/sub-025/ses-meg01/sub25ses01_3011085.02_20170118_01.ds';
subjects(25).logfile            = '/project/3011085.02/raw/sub-025/ses-beh01/sub25ses01.mat';
subjects(25).headshape          = '/project/3011085.02/raw/sub-025/polhemus_erfosc_sub-025.pos';
subjects(25).mri                = '/project/3011085.02/raw/sub-025/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.43.67027.2017011814031851926362919.IMA';
subjects(25).mridir             = '/project/3011085.02/processed/sub-025/ses-mri01/';
subjects(25).ecgcomp            = [2, 6, 14];
subjects(25).eyecomp            = [1];
subjects(25).badtrials          = [];
subjects(25).gammaband          = [38 58];%TFR [34 60];powspec
subjects(25).gammapeak          = 48;

subjects(26).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(26).dataset            = '/project/3011085.02/raw/sub-026/ses-meg01/sub26ses01_3011085.02_20170118_01.ds';
subjects(26).logfile            = '/project/3011085.02/raw/sub-026/ses-beh01/sub26ses01.mat';
subjects(26).headshape          = '/project/3011085.02/raw/sub-026/polhemus_erfosc_sub-026.pos';
subjects(26).mri                = '/project/3011085.02/raw/sub-026/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/PAUGAA_20161205_S15.MR.DCCN_SKYRA.0002.0001.2016.12.05.10.25.31.305657.1318410308.IMA';
subjects(26).mridir             = '/project/3011085.02/processed/sub-026/ses-mri01/';
subjects(26).ecgcomp            = [2, 6, 11];
subjects(26).eyecomp            = [1];
subjects(26).badtrials          = [1, 399];
subjects(26).gammaband          = [42 54]%TFR[42 52];powspec
subjects(26).gammapeak          = 48;

subjects(27).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(27).dataset            = '/project/3011085.02/raw/sub-027/ses-meg01/sub027ses01_3011085.02_20170123_01.ds';
subjects(27).logfile            = '/project/3011085.02/raw/sub-027/ses-beh01/sub27ses01.mat';
subjects(27).headshape          = '/project/3011085.02/raw/sub-027/polhemus_erfosc_sub-027.pos';
subjects(27).mri                = '/project/3011085.02/raw/sub-027/ses-mri01/002-t1_mpr_ns_PH8/PAUGAA_20140923_S20.MR.DCCN_SEQUENCES_STANDARD_SEQUENCES.0002.0001.2014.09.23.14.52.22.958467.43979878.IMA';
subjects(27).mridir             = '/project/3011085.02/processed/sub-027/ses-mri01/';
subjects(27).ecgcomp            = [1, 22];
subjects(27).eyecomp            = [5];
subjects(27).badtrials          = [];
subjects(27).gammaband          = [58 70];%TFR [50 78];powspec
subjects(27).gammapeak          = 66;

subjects(28).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(28).dataset            = '/project/3011085.02/raw/sub-028/ses-meg01/sub028ses01_3011085.02_20170125_01.ds';
subjects(28).logfile            = '/project/3011085.02/raw/sub-028/ses-beh01/sub28ses01.mat';
subjects(28).headshape          = '/project/3011085.02/raw/sub-028/polhemus_erfosc_sub-028.pos';
subjects(28).mri                = '/project/3011085.02/raw/sub-028/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.19.45416.2017012510404133866953306.IMA';
subjects(28).mridir             = '/project/3011085.02/processed/sub-028/ses-mri01/';
subjects(28).ecgcomp            = [5, 21];
subjects(28).eyecomp            = [1];
subjects(28).badtrials          = [190, 453];
subjects(28).gammaband          = [34 46];%TFR[30 44; 68 86]; %powspec
subjects(28).gammapeak          = [34; 76]; % 34 strongest. 

subjects(29).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(29).dataset            = '/project/3011085.02/raw/sub-029/ses-meg01/sub029ses01_3011085.02_20170303_01.ds';
subjects(29).logfile            = '/project/3011085.02/raw/sub-029/ses-beh01/sub29ses01.mat';
subjects(29).headshape          = '/project/3011085.02/raw/sub-029/polhemus_erfosc_sub-029.pos';
subjects(29).mri                = '/project/3011085.02/raw/sub-029/ses-mri01/015-t1_mprage_sag_p2_iso_1.0/BRAZAN_SUBJ28_20161006.MR.BRAZAN_SKYRA.0015.0001.2016.10.06.15.51.35.148367.1280435354.IMA';
subjects(29).mridir             = '/project/3011085.02/processed/sub-029/ses-mri01/';
subjects(29).ecgcomp            = [3, 15];
subjects(29).eyecomp            = [1, 2];
subjects(29).badtrials          = [];
subjects(29).gammaband          = [50 66];%TFR [40 72];%powspec
subjects(29).gammapeak          = 58; 

subjects(30).channels           = {'MEG', '-MRT32', 'HLC*', '-MRC11', '-MRF61'};
subjects(30).dataset            = '/project/3011085.02/raw/sub-030/ses-meg01/sub030ses01_3011085.02_20170125_01.ds';
subjects(30).logfile            = '/project/3011085.02/raw/sub-030/ses-beh01/sub30ses01.mat';
subjects(30).headshape          = '/project/3011085.02/raw/sub-030/polhemus_erfosc_sub-030.pos';
subjects(30).mri                = '/project/3011085.02/raw/sub-030/ses-mri01/002-t1_mprage_sag_p2_iso_1.0/20170105_CHRUTZ_MEGS02.MR.CHRUTZ_PRISMA.0002.0001.2017.01.05.13.46.41.183639.604120162.IMA';
subjects(30).mridir             = '/project/3011085.02/processed/sub-030/ses-mri01/';
subjects(30).ecgcomp            = [3, 9, 14];
subjects(30).eyecomp            = [1, 2];
subjects(30).badtrials          = []; 
subjects(30).gammaband          = [50 68];%TFR [48 88];powspec %absent? also looks like it in TFR
subjects(30).gammapeak          = 74;

subjects(31).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(31).dataset            = '/project/3011085.02/raw/sub-031/ses-meg01/sub031ses01_3011085.02_20170126_01.ds';
subjects(31).logfile            = '/project/3011085.02/raw/sub-031/ses-beh01/sub31ses01.mat';
subjects(31).headshape          = '/project/3011085.02/raw/sub-031/polhemus_erfosc_sub-031.pos';
subjects(31).mri                = '/project/3011085.02/raw/sub-031/ses-mri01/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.43.66068.2017012710350118315501039.IMA';
subjects(31).mridir             = '/project/3011085.02/processed/sub-031/ses-mri01/';
subjects(31).ecgcomp            = [3, 7];
subjects(31).eyecomp            = [2, 10];
subjects(31).badtrials          = [];
subjects(31).gammaband          = [50 66];%[46 70];
subjects(31).gammapeak          = 58;

subjects(32).channels           = {'MEG', '-MRT32', 'HLC*'};
subjects(32).dataset            = '/project/3011085.02/raw/sub-032/ses-meg01/sub031ses01_3011085.02_20170126_02.ds';
subjects(32).logfile            = '/project/3011085.02/raw/sub-032/ses-beh01/sub32ses01.mat';
subjects(32).headshape          = '/project/3011085.02/raw/sub-032/polhemus_erfosc_sub-032.pos';
subjects(32).mri                = '/project/3011085.02/raw/sub-032/ses-mri02/002-t1_mprage_sag_p2_iso_1.0_20ch_head/00001_1.3.12.2.1107.5.2.43.66068.2017021514302695199429043.IMA';
subjects(32).mridir             = '/project/3011085.02/processed/sub-032/ses-mri02/';
subjects(32).ecgcomp            = [3, 20, 24];
subjects(32).eyecomp            = [1, 13];
subjects(32).badtrials          = [];
subjects(32).gammaband          = [46 58];%TFR %[34 64];powspec
subjects(32).gammapeak          = 52;

subjects(33).channels           = {'MEG', '-MRT32', '-MLF12', 'HLC*'};
subjects(33).dataset            = '/project/3011085.02/raw/sub-033/ses-meg01/sub033ses01_3011085.02_20170126_01.ds';
subjects(33).logfile            = '/project/3011085.02/raw/sub-033/ses-beh01/sub33ses01.mat';
subjects(33).headshape          = '/project/3011085.02/raw/sub-033/polhemus_erfosc_sub-033.pos';
subjects(33).mri                = '/project/3011085.02/raw/sub-033/ses-mri01/20160612_125811t1mpragesagp2iso10s010a1001.nii.gz';
subjects(33).mridir             = '/project/3011085.02/processed/sub-033/ses-mri01/';
subjects(33).ecgcomp            = [2, 6];
subjects(33).eyecomp            = [1];
subjects(33).badtrials          = [12, 27, 181, 207];
subjects(33).gammaband          = [54 68];%powspec + TFR
subjects(33).gammapeak          = 60;
