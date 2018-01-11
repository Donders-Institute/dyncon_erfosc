function [source_parc] = erfosc_lcmv_parc(data_shift, headmodel, sourcemodel)

cfg         = [];
cfg.latency = [-0.1 inf];
data_shift  = ft_selectdata(cfg, data_shift);
cfg.latency = [-0.05 0.25];
data_shift_short = ft_selectdata(cfg, data_shift);

cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.1 0];
cfg.covariance = 'yes';
cfg.removemean = 'no';
tlck = ft_timelockanalysis(cfg, data_shift_short);
tlck_long = ft_timelockanalysis(cfg, data_shift);

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

tmp     = rmfield(data, {'elec' 'grad'});
selparc = setdiff(1:numel(atlas.parcellationlabel),[1 2 194 195]); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck.time;
source_parc.F     = cell(numel(source_parc.label),1);
source_parc.avg   = zeros(numel(selparc),numel(source_parc.time));
source_parc.dimord = 'chan_time';

for k = 1:numel(selparc)
  tmpF = F(atlas.parcellation==selparc(k),:);
  tmp.trial = tmpF*data.trial;
  tmp.label = data.label(1:size(tmpF,1));
  tmpcomp   = ft_componentanalysis(cfg, tmp);

  source_parc.F{k}     = tmpcomp.unmixing*tmpF;
  source_parc.avg(k,:) = source_parc.F{k}(1,:)*tlck.avg;
end





