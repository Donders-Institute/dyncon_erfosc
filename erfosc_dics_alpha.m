function [source_onset, source_shift, Tval, F] = erfosc_dics_alpha(freq_onset, freq_shift, headmodel, sourcemodel)
% Compute power in the low frequencies on the source level using a DICS
% beamformer (Gross et al, 2001). Also estimates the location with the 
% highest increase from baseline.
%
% INPUT 
%   freq_onset: frequency estimate, timing relative to stimulus onset
%   freq_shift: frequency estimate, timing relative to stimulus change
%   headmodel: volume conduction model of the head, output of
%       ft_prepare_headmodel
%   sourcemodel: 2D or 3D sourcemodel
% 
% OUTPUT
%   source_onset: source estimate, timing relative to stimulus onset
%   source_shift: source estimate, timing relative to stimulus change
%   Tval: T values in source space of the comparisson vs baseline
%   F: spatial filters
cfg           = [];
cfg.appenddim = 'rpt';
cfg.parameter = 'fourierspctrm';
freq          = ft_appendfreq([], freq_onset, freq_shift);

if ~isfield(sourcemodel, 'leadfield')
  cfg           = [];
  cfg.headmodel = ft_convert_units(headmodel, 'm');
  cfg.grid      = ft_convert_units(sourcemodel, 'm');
  cfg.grad      = ft_convert_units(freq_onset.grad, 'm');
  cfg.channel   = freq_onset.label;
  sourcemodel   = ft_prepare_leadfield(cfg);
end

cfg = [];
cfg.method = 'dics';
cfg.frequency = freq_onset.freq(1);
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.dics.keepfilter = 'yes';
cfg.dics.fixedori   = 'yes';
cfg.dics.realfilter = 'yes';
cfg.dics.lambda     = '100%';
source = ft_sourceanalysis(cfg, ft_checkdata(freq, 'cmbrepresentation', 'fullfast'));
cfg.grid.filter = source.avg.filter;
cfg.dics.keepfilter = 'no';

source_onset = ft_sourceanalysis(cfg, ft_checkdata(freq_onset, 'cmbrepresentation', 'fullfast'));
source_shift = ft_sourceanalysis(cfg, ft_checkdata(freq_shift, 'cmbrepresentation', 'fullfast'));

% projection matrix to get from fourier to power
nrpt = numel(freq_onset.cumtapcnt);
ntap = freq_onset.cumtapcnt(1);

ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
iy = 1:(nrpt*ntap);
iz = ones(nrpt*ntap,1)./ntap;
P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);

% hacky T-statistic computation
F = cat(1, source.avg.filter{:});
p_onset = (abs(F*transpose(freq_onset.fourierspctrm)).^2)*P;
p_shift = (abs(F*transpose(freq_shift.fourierspctrm)).^2)*P;

cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
design = [ones(1,nrpt)*2 ones(1,nrpt);1:nrpt 1:nrpt];
Tstat_dep = ft_statfun_depsamplesT(cfg, [p_onset p_shift], design);
cfg = [];
cfg.ivar = 1;
design = design(1,:);
Tstat_indep = ft_statfun_indepsamplesT(cfg, [p_onset p_shift], design);

Tval = zeros(size(source.pos,1),1);
Tval(source.inside) = Tstat_dep.stat;

Ftmp = zeros(size(source_onset.pos,1), size(F,2));
Ftmp(source_onset.inside, :) = F;
F = Ftmp; clear Ftmp;

