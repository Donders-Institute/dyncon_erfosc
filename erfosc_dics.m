function [sourcefreq_onset, sourcefreq_shift, Tval] = erfosc_dics(freq_onset, freq_shift, headmodel, sourcemodel)
% Compute power on the source level using a DICS
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
  headmodel      = ft_convert_units(headmodel, 'm');
  sourcemodel    = ft_convert_units(sourcemodel, 'm');

  cfg = [];
  cfg.comment = 'convert to unit `m`.';
  headmodel   = ft_annotate(cfg, headmodel);
  sourcemodel = ft_annotate(cfg, sourcemodel);
  
  cfg           = [];
  cfg.headmodel = headmodel;
  cfg.grid      = sourcemodel;
  cfg.grad      = ft_convert_units(freq_onset.grad, 'm');
  cfg.channel   = freq_onset.label;
  cfg.comment   = 'get grad and label information from frequency data.';
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
cfg.comment = 'change frequency data to cmbrepresentation `fullfast` with ft_checkdata.';
source = ft_sourceanalysis(cfg, ft_checkdata(freq, 'cmbrepresentation', 'fullfast'));
cfg.grid.filter = source.avg.filter;

cfg.comment = 'change frequency data to cmbrepresentation `fullfast` with ft_checkdata. Take the spatial filter from the previous step.';
sourcefreq_onset = ft_sourceanalysis(cfg, ft_checkdata(freq_onset, 'cmbrepresentation', 'fullfast'));
sourcefreq_shift = ft_sourceanalysis(cfg, ft_checkdata(freq_shift, 'cmbrepresentation', 'fullfast'));

% projection matrix to get from fourier to power
nrpt = numel(freq_onset.cumtapcnt);
ntap = freq_onset.cumtapcnt(1);

ix = reshape(repmat(1:nrpt,[ntap 1]),[],1);
iy = 1:(nrpt*ntap);
iz = ones(nrpt*ntap,1)./ntap;
P  = sparse(iy,ix,iz,nrpt*ntap,nrpt);

% T-statistic computation
F = zeros(numel(source.inside), numel(freq_onset.label));
F(source.inside,:) = cat(1, source.avg.filter{:}); %FIXME do I need to change this into a FT structure?

cfg=[];
cfg.comment = 'concatenate the spatial filters';
sourcefreq_onset.F = F;
sourcefreq_onset = ft_annotate(cfg, sourcefreq_onset);
sourcefreq_shift.F = F;
sourcefreq_shift = ft_annotate(cfg, sourcefreq_shift);

p_onset = rmfield(source, 'avg');
p_onset.dimord = 'pos_rpt';
p_onset.pow = (abs(F*transpose(freq_onset.fourierspctrm)).^2)*P;
cfg=[];
cfg.comment = 'calculate single-trial source power (locked to stimulus onset) by computing "pow = (abs(spatial_filter*transpose(fourierspctrm)).^2)*P", where P is the weigth matrix for combining over tapers.';
p_onset = ft_annotate(cfg, p_onset);

p_shift = rmfield(source, 'avg');
p_shift.dimord = 'pos_rpt';
p_shift.pow = (abs(F*transpose(freq_shift.fourierspctrm)).^2)*P;
cfg.comment = 'calculate single-trial source power (locked to stimulus change) by computing "pow = (abs(spatial_filter*transpose(fourierspctrm)).^2)*P", where P is the weigth matrix for combining over tapers.';
p_shift = ft_annotate(cfg, p_shift);

cfg = [];
cfg.ivar = 1;
cfg.uvar = 2;
design = [ones(1,nrpt)*2 ones(1,nrpt);1:nrpt 1:nrpt];
cfg.design = design;
cfg.method = 'analytic';
cfg.statistic = 'depsamplesT';
Tval = ft_sourcestatistics(cfg, p_onset, p_shift);


