function [source_parc] = erfosc_lcmv_parc(data_shift, headmodel, sourcemodel, doresplocked)
load('atlas_subparc374_8k.mat')
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
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
cfg.lcmv.lambda     = '100%';
cfg.lcmv.weightnorm = 'unitnoisegain';
source = ft_sourceanalysis(cfg, tlck);
F      = zeros(size(source.pos,1),numel(tlck.label));
F(source.inside,:) = cat(1,source.avg.filter{:});

% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

tmp     = rmfield(data_shift, {'elec' 'grad'});
exclude_label = match_str(atlas.parcellationlabel, {'L_???_01', 'L_MEDIAL.WALL_01', 'R_???_01', 'R_MEDIAL.WALL_01'});
selparc = setdiff(1:numel(atlas.parcellationlabel),exclude_label); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck.time;
source_parc.F     = cell(numel(source_parc.label),1);
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

  source_parc.F{k}     = tmpcomp.unmixing*tmpF;
  source_parc.avg(k,:) = source_parc.F{k}(1,:)*tlck.avg;
end





