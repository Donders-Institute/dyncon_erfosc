function erf_osc_analysis_tfa_allsubs(freqRange, zeropoint)
% Loads analysis from erf_osc_analysis_tfa, makes a relative contrast with
% baseline and averages over subjects.
%
% INPUT
%   freRange (string): 'low' (default) or 'high', frequency range (below or 
%       above 30 Hz; affects configuration settings).
%   zeropoint (string): 'onset' (default) or 'reversal', what to time lock 
%       the data to. if 'onset' a seperate baseline estimate is saved.
% 
% OUTPUT
%   data saved on disk
%   d: group-average relative difference between stimulus and baseline
%       conditions.

if ~exist(freqRange); freqRange = 'low'; end
if ~exist(zeropoint); zeropoint = 'onset'; end

erf_osc_datainfo;
k=1;
for subj=allsubs
    tmp = load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_tfa_%s_%s', subj, subj, freqRange, zeropoint), 'tfa');
    tfa{k} = tmp.tfa;
    tmp = load(sprintf('/project/3011085.02/analysis/freq/sub-%03d/sub-%03d_tfa_%s_%s', subj, subj, freqRange, 'onset'), 'baseline');
    baseline{k} = tmp.baseline;
    k=k+1;
end

for k=1:numel(tfa);
    baseline{k}.powspctrm = repmat(mean(baseline{k}.powspctrm,3),[1 1 numel(tfa{k}.time)]);
    tfa{k}.powspctrm = tfa{k}.powspctrm./baseline{k}.powspctrm;
end

cfg=[];
cfg.appenddim = 'rpt';
d = ft_appendfreq(cfg, tfa{:});

save(sprintf('/project/3011085.02/analysis/tfr_all_%s_%s.mat', freqRange, zeropoint), 'd', '-v7.3');



