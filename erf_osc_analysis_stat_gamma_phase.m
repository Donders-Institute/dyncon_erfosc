option = 'tinput';
% option = 20;

erf_osc_datainfo;

if ~ischar(option)
    option=num2str(option);
end
for subj=allsubs
tmp{subj} = load(sprintf('/project/3011085.02/results/freq/sub-%03d/gamma_angle_%s.mat', subj, option));
end
k=1;
for subj=allsubs
stat{subj} = tmp{subj}.stat;
statrand{subj} = mean(tmp{subj}.statrand, 2);
inputTime(k) = tmp{subj}.inputTime;
k=k+1;
if subj==2;
dum = tmp{subj}.statrand;
end
end
stat_all = (cat(2, stat{allsubs}))';
statrand_all = (cat(2, statrand{allsubs}))';

t=0:1/1200:0.5;

for k=1:32
idxT(k) = nearest(t, inputTime(k));
end

for k=1:32
d(k) = stat_all(k, idxT(k))./mean(statrand_all(k, :),2);
end