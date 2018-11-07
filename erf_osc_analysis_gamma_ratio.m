% estimates the gamma power increase (ratio) from baseline at the
% individual gamma peak frequency. This estimation is based on either the
% average over occipital channels, or the maximum occipital channel.

erf_osc_datainfo;

k=1;
for subj=allsubs
tmp{k} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/pow.mat', subj));
k=k+1;
end


for k=1:32
peakFreq(k) = tmp{k}.peakFreq_gamma;
ratio_broadband(k) = tmp{k}.gamRatio; % mean power increase over channel and frequency
powRatio{k} = tmp{k}.powRatio;
end


for k=1:32
f(k) = find(powRatio{k}.freq==peakFreq(k));
ratio_maxchan(k) = max(powRatio{k}.powspctrm(:,f(k))); % power at max channel at peak frequency
end

ratio_maxchan_GA = mean(ratio_maxchan);

%% reaction times
k=1;
for subj=allsubs
    load(sprintf('/project/3011085.02/results/behavior/sub-%03d/rt.mat', subj));
    RT{k} = rt;
    clear rt;
    k=k+1;
end

%% correlations with age.
load('/project/3011085.02/age.mat');