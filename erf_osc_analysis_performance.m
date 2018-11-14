function erf_osc_analysis_perfomance(subj)

load(sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_trialinfo.mat', subj, subj));
trialinfo_all = trialinfo;
nTrials_all = size(trialinfo,1);
idxNoReversal = find(trialinfo(:,8)==0);
noReversal = length(idxNoReversal);
trialinfo(idxNoReversal,:)=[];
nTrials_reversal = size(trialinfo,1);
nTrials_validResp = sum(trialinfo(:,9)>0 & ((trialinfo(:,9)-trialinfo(:,8))/1200)<0.7);
performanceAll = nTrials_validResp/nTrials_reversal;
clear trialinfo


load(sprintf('/project/3011085.02/processed/sub-%03d/ses-meg01/sub-%03d_cleandata.mat', subj, subj), 'dataClean');
trialinfo_cleaned = dataClean.trialinfo;
trialinfo = trialinfo_cleaned;
nTrials_cleaned = size(trialinfo,1);
idxNoReversal_cleaned = find(trialinfo(:,5)==0);
noReversal_cleaned = length(idxNoReversal_cleaned);
trialinfo(idxNoReversal_cleaned,:)=[];
nTrials_reversal_cleaned = size(trialinfo,1);
nTrials_validResp_cleaned = sum(trialinfo(:,6)>0 & ((trialinfo(:,6)-trialinfo(:,5))/1200)<0.7);
performance_cleaned = nTrials_validResp_cleaned/nTrials_reversal_cleaned;

save(sprintf('/project/3011085.02/analysis/behavior/sub-%03d/sub-%03d_performance.mat', subj, subj), ...
    'performanceAll', 'nTrials_validResp', 'nTrials_all', 'nTrials_reversal', 'noReversal', 'trialinfo_all', ...
    'trialinfo_cleaned', 'nTrials_cleaned', 'idxNoReversal_cleaned', 'noReversal_cleaned', 'nTrials_reversal_cleaned', ...
    'nTrials_validResp_cleaned', 'performance_cleaned');
end


