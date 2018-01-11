function [freq_onset, freq_shift, P] = erfosc_freq(data_onset, data_shift, latency, subject, foi, smo)

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
