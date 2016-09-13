% function headlocation(subj)

% define trials

datainfo;  
for subj=1:10
    jj=1;
for ses = subjects(subj).sessions
cfg = [];
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);
 
% preprocess the headposition data
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                              'HLC0021','HLC0022','HLC0023', ...
                              'HLC0031','HLC0032','HLC0033'};
cfg.continuous = 'yes'; 

headpos = ft_preprocessing(cfg);
 
% calculate the mean coil position per trial
ntrials = length(headpos.sampleinfo);
for t = 1:ntrials
coil1{subj}.ses{jj}(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
coil2{subj}.ses{jj}(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
coil3{subj}.ses{jj}(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
end
jj=jj+1;
end
end
%
% addpath(genpath('/home/common/matlab/fieldtrip/test/private/circumcenter.m')); % make sure the multivariate toolbox is in your path
% calculate the headposition and orientation per trial
for subj = 1:10
    jj=1;
    kk=1;
    for ses = subjects(subj).sessions
        ntrial_ses = length(coil2{subj}.ses{kk});
    coil1_sessions(:,jj:jj+ntrial_ses-1) = coil1{subj}.ses{kk};
    coil2_sessions(:,jj:jj+ntrial_ses-1) = coil2{subj}.ses{kk};
    coil3_sessions(:,jj:jj+ntrial_ses-1) = coil3{subj}.ses{kk};
    jj = jj+ntrial_ses;
    kk=kk+1;
    end

cc{subj} = circumcenter(coil1_sessions, coil2_sessions, coil3_sessions);
cc_mean{subj} = mean(cc{subj},2); 
% demean to obtain translations and rotations from the average position and orientation
cc_dem{subj} = [cc{subj} - repmat(cc_mean{subj},1,size(cc{subj},2))]';
end
jj=1;
for subj=1:10
avg_cc_dem(:,jj) = mean(cc_dem{subj}',2);
jj=jj+1;
end


for jj=1:10
vector(jj) = sqrt(avg_cc_dem(1,jj)^2+avg_cc_dem(2,jj)^2+avg_cc_dem(3,jj)^2);
end
bar(vector)
figure; hist(vector,10)

%%

%{
datainfo;
sessions = subjects(subj).sessions;
allcc_rel=[];
length_ses = zeros(1, length(sessions));
for ii=1:length(sessions)
  ses = sessions{ii};
  
cfg                         = [];
cfg = subjects(subj);
cfg.dataset = subjects(subj).session(ses).dataset;
cfg.trialfun = subjects(subj).trialfun;
cfg = ft_definetrial(cfg);


cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
                              'HLC0021','HLC0022','HLC0023', ...
                              'HLC0031','HLC0032','HLC0033'};
cfg.continuous = 'yes'; 
headpos = ft_preprocessing(cfg);
% coil calculation
% calculate the mean coil position per trial
ntrials = length(headpos.sampleinfo);
for t = 1:ntrials
coil1(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
coil2(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
coil3(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
end

cc = ft_circumcenter(coil1, coil2, coil3);
cc_dem = [cc - repmat(mean(cc,2),1,size(cc,2))]';
if ii==1
  length_ses(ii) = length(cc);
else
  length_ses(ii) = length(cc)+length_ses(ii-1);
end

cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))];
icc_rel=cc_rel';
maxposchange = max(abs(icc_rel(:,1:3)*1000)) % in mm 
allcc_rel = [allcc_rel cc_rel];
end 
allcc_rel=allcc_rel';
% plot translations
figure(); 
hold on
plot(allcc_rel(:,1:3)*1000) % in mm
for ii=1:length(sessions)
  vline(length_ses(ii), 'k')
end
legend('x', 'y','z')
title('head translation relative to first trial of first session')
xlabel('trial')
ylabel('movement (mm)')
 
% plot rotations
figure(); 
hold on
plot(allcc_rel(:,4:6))
for ii=1:length(sessions)
  vline(length_ses(ii), 'k')
end

title('head rotation relative to first trial of first session')
xlabel('trial')
ylabel('angle (degrees)')
legend('x', 'y','z')
% end
%}