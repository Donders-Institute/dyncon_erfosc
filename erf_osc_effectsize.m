
%% GLM gamma
load('/project/3011085.02/results/stat_glm_gamma_time_reversal.mat')
x1 = betas_plcmb_GA.trial.*(permute(repmat(stat.posclusterslabelmat==1,[1 1 32]),[3,1,2]));
x2 = betas_bl_plcmb_GA.trial.*(permute(repmat(stat.posclusterslabelmat==1,[1 1 32]),[3,1,2]));
x1 = (x1(find(x1(:)~=0)));
x2 = (x2(find(x2(:)~=0)));
paired_sd = std(x1-x2);
cohensd = (mean(x1-x2))/paired_sd


%% or???

x1 = betas_plcmb_GA.trial;
x2 = betas_bl_plcmb_GA.trial;
x1 = x1(:,stat.posclusterslabelmat==1);
x2 = x2(:,stat.posclusterslabelmat==1);
X1 = mean(x1,2);
X2 = mean(x2,2);
d = (mean(X1-X2))/(std(X1-X2))
