function [freq_onset, freq_shift] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo)
% Frequency analysis using fieldtrip mtmfft
%
% INPUT
%   data_onset (struct): fieldtrip data structure, time locked to stimulus
%       onset.
%   data_shift (struct): fieldtrip data structure, time locked to stimulus
%       change.
%   latency: scalar or string, can be 'all', 'prestim', 'poststim', or 
%       [beg end], specify time range in seconds.
%   subject (int): subject ID, ranging from 1 to 33, excluding 10.
%   foi: vector 1 x numfoi, frequencies of interest (default: subject-
%       specific gamma peak frequency)
%   smo: number, the amount of spectral smoothing through multi-tapering. 
%       Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz 
%       smoothing box. (default: subject specific gamma bandwidth)
% 
% OUTPUT
%   freq_onset: frequency estimate, timing relative to stimulus onset
%   freq_shift: frequency estimate, timing relative to stimulus change
%   P: projection matrix to go from fourier to power


if nargin<5 || isempty(foi)
  foi = [1 1].*subject.gammapeak(end);
end
if nargin<6 || isempty(smo)
  smo = diff(subject.gammaband(end,:))./2;
end

cfg         = [];
cfg.latency = latency;
data_onset  = ft_selectdata(cfg, data_onset);
data_shift  = ft_selectdata(cfg, data_shift);

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = foi;
cfg.tapsmofrq = smo;
cfg.pad    = 1;
freq_onset = ft_freqanalysis(cfg, data_onset);
freq_shift = ft_freqanalysis(cfg, data_shift);

% projection matrix to get from fourier to power
nrpt = numel(freq_onset.cumtapcnt);
ntap = freq_onset.cumtapcnt(1);

ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
iy = 1:(nrpt*ntap);
iz = ones(nrpt*ntap,1)./ntap;
P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);

freq_onset.fourier2pow = P;
freq_shift.fourier2pow = P;
cfg=[];
cfg.comment = 'add projection matrix (field P) to get from fourier to power';
freq_onset = ft_annotate(cfg, data_onset);
freq_shift = ft_annotate(cfg, data_shift);
