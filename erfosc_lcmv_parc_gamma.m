function [source_parc] = erfosc_lcmv_parc_gamma(data_shift, headmodel, sourcemodel, atlas, foi)

cfg=[];
cfg.bpfilter = 'yes';
cfg.bpfreq = foi;
cfg.bpfilttype = 'firws';
data_shift_bp = ft_preprocessing(cfg, data_shift);

cfg         = [];
cfg.latency = [-0.75 0];
data_pre_short  = ft_selectdata(cfg, data_shift_bp);
cfg.latency = [-1.5 0.5];
data_pre = ft_selectdata(cfg, data_shift);

cfg = [];
cfg.covariance = 'yes';
cfg.removemean = 'yes';
tlck_pre_short = ft_timelockanalysis(cfg, data_pre_short);
tlck_pre = ft_timelockanalysis(cfg, data_pre);

if ~isfield(sourcemodel, 'leadfield')
  cfg           = [];
  cfg.headmodel = ft_convert_units(headmodel, 'm');
  cfg.grid      = ft_convert_units(sourcemodel, 'm');
  cfg.grad      = ft_convert_units(data_pre_short.grad, 'm');
  cfg.channel   = data_pre_short.label;
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
source = ft_sourceanalysis(cfg, tlck_pre_short); % spatial filter based on pre shift
F      = zeros(size(source.pos,1),numel(tlck_pre_short.label));
F(source.inside,:) = cat(1,source.avg.filter{:});

% prepare the cfg for pca
cfg                       = [];
cfg.method                = 'pca';

tmp     = rmfield(data_pre_short, {'elec' 'grad'});
selparc = setdiff(1:numel(atlas.parcellationlabel),[1 2 194 195]); % hard coded exclusion of midline and ???

source_parc.label = atlas.parcellationlabel(selparc);
source_parc.time  = tlck_pre.time;
source_parc.F     = cell(numel(source_parc.label),1);
source_parc.avg   = zeros(numel(selparc),numel(source_parc.time));
source_parc.dimord = 'chan_time';

for k = 1:numel(selparc)
  tmpF = F(atlas.parcellation==selparc(k),:);
  tmp.trial = tmpF*data_pre_short.trial; % component analysis based on short pre shift window
  for l=1:size(tmpF,1)
      tmplabel{l} = sprintf('loc%03d', l);
  end
  tmp.label = tmplabel;
  clear tmplabel
  tmpcomp   = ft_componentanalysis(cfg, tmp);

  source_parc.F{k}     = tmpcomp.unmixing*tmpF;
end





