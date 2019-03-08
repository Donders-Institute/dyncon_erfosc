function [source_parc] = erfosc_lcmv_parc(data_shift, headmodel, sourcemodel, atlas, doresplocked)
% Model time courses on the source level using a LCMV beamformer. By
% default the evoked activity is estimated at stimulus change, but can be
% done at response onset. dipoles are combined into anatomically defined
% parcels.
%
% INPUT 
%   data_shift: data, timing relative to stimulus change
%   headmodel: volume conduction model of the head, output of
%       ft_prepare_headmodel
%   sourcemodel: 2D or 3D sourcemodel
%   atlas: anatomical or functional atlas describing how to combine dipole
%      locations. Only implemented for HBM anatomical atlas (374parc)
%   doresplocked (logical): true when lcmv should be optimized for response
%       locked activity.
% 
% OUTPUT
%   source_parc: channel level fieldtrip data structure, containing single
%       trial time courses for each parcel defined in the atlas.


if ~exist('doresplocked'); doresplocked=false; end

cfg         = [];
if doresplocked
    cfg.latency = [-0.5 0.4];
else
    cfg.latency = [-0.1 0.6];
end
data_shift = ft_selectdata(cfg, data_shift);

cfg = [];
cfg.preproc.demean = 'yes';
if ~doresplocked
    cfg.preproc.baselinewindow = [-0.1 0];
    cfg.removemean = 'no';
end
cfg.covariance = 'yes';

tlck = ft_timelockanalysis(cfg, data_shift);

if ~isfield(sourcemodel, 'leadfield')
  cfg           = [];
  cfg.headmodel = ft_convert_units(headmodel, 'm');
  cfg.grid      = ft_convert_units(sourcemodel, 'm');
  cfg.grad      = ft_convert_units(data_shift.grad, 'm');
  cfg.channel   = data_shift.label;
  sourcemodel   = ft_prepare_leadfield(cfg);
end

cfg = [];
cfg.method = 'lcmv';
cfg.headmodel = headmodel;
cfg.grid      = sourcemodel;
cfg.keepleadfield = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '100%';
%cfg.lcmv.weightnorm = 'unitnoisegain'; % this confuses me in terms of the
%unit-gain inspection when keeping the leadfields: it's also just a scaling
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepori = 'yes';
source = ft_sourceanalysis(cfg, tlck);
F      = zeros(size(source.pos,1),numel(tlck.label));
F(source.inside,:) = cat(1,source.avg.filter{:});

L      = zeros(numel(tlck.label), size(source.pos,1));
L(:,source.inside) = cat(2,source.leadfield{:});

% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

tmp     = rmfield(data_shift, {'elec' 'grad'});
selparc = 1:numel(atlas.parcellationlabel);
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck.time;
source_parc.F     = cell(numel(source_parc.label),1);
source_parc.L     = cell(numel(source_parc.label),1);
source_parc.avg   = zeros(numel(selparc),numel(source_parc.time));
source_parc.dimord = 'chan_time';

for k = 1:numel(selparc)
  tmpF = F(atlas.parcellation==selparc(k),:);
  tmp.trial = tmpF*data_shift.trial;
  for l=1:size(tmpF,1)
      tmplabel{l} = sprintf('loc%03d', l);
  end
  tmp.label = tmplabel;
  clear tmplabel
  tmpcomp   = ft_componentanalysis(cfg, tmp);

  tmpL = L(:,atlas.parcellation==selparc(k));
  
  source_parc.F{k}     = tmpcomp.unmixing*tmpF;
  source_parc.L{k}     = tmpL*tmpcomp.unmixing';
  source_parc.avg(k,:) = source_parc.F{k}(1,:)*tlck.avg;
end





