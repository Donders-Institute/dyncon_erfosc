function erf_osc_analysis_corr(subj, isPilot)

if nargin<1 || isempty(subj)
    subj = 1;
end
if nargin<2 || isempty(isPilot)
    isPilot = false;
end


% initiate diary
ft_diary('on')
erf_osc_datainfo;

%% gamma power - reaction time

for subj=allsubs
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_peak.mat', subj), 'peakFreq_gamma');
    load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_virtual_channel.mat', subj), 'gammaChan'); % load gamma power
    rt{subj} = load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj)); % load gamma power
    
    for i=1:length(gammaChan.trial)
        gammaPow_tmp(i) = log(gammaChan.trial(i).pow);
    end
    gammaPow{subj} = gammaPow_tmp-mean(gammaPow_tmp);
    clear gammaChan gammaPow_tmp
    rt{subj} = rt{subj}.rt;
    [r(subj) p(subj)] = corr(gammaPow{subj}', rt{subj}, 'type', 'spearman');
end


%% P1 amplitude - TFR



%% save
if strcmp(correlation, 'gamma_rt');
    if isPilot
        filename = sprintf('/project/3011085.02/results/freq/pilot-%03d/corr_gamma_rt', subj);
    else
        filename = sprintf('/project/3011085.02/results/freq/sub-%03d/corr_gamma_rt', subj);
    end
    save(fullfile([filename '.mat']), 'r', 'p', 'rt', 'gammaPow');
end
ft_diary('off')
end
