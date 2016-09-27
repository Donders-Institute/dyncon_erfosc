%% Plot topography

cfg=[];
cfg.layout = 'CTF275_helmet.mat';
cfg.latency = [0 0.2];
cfg.xlim=[0 0.2];
cfg.zlim=[-1e-13 1e-13];
cfg.colormap = 'jet';

% based on grating onset (baseline window corrected)
figure;
subplot(1,3,1);ft_topoplotER(cfg, tlOnsetM_rs);title('Middle, onset (baseline window corrected)');
% based on grating shift (baseline window corrected)
subplot(1,3,2);ft_topoplotER(cfg, tlShiftM_rs);title('Middle, shift (baseline window corrected)');
% based on grating shift (100ms preShift corrected)
subplot(1,3,3);ft_topoplotER(cfg, tlShiftPreM_rs);title('Middle, shift (baseline 100ms pre)');

%% plot timecourses

cfg=[];
cfg.xlim = [-0.5 1];
cfg.ylim = [-5e-13 3e-13];
figure;
% based on grating onset (baseline window corrected)
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlOnsetL_rs);title('Left, onset (baseline window corrected)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlOnsetM_rs);title('Middle, onset (baseline window corrected)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlOnsetR_rs);title('Right, onset (baseline window corrected)');
% based on grating shift (baseline window corrected)
figure;
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlShiftL_rs);title('Left, shift (baseline window corrected)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlShiftM_rs);title('Middle, shift (baseline window corrected)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlShiftR_rs);title('Right, shift (baseline window corrected)');
% based on grating shift (200ms preShift corrected)
figure;
cfg.channel = 'MRO11'; % negative peak in topo left
subplot(1,3,1);ft_singleplotER(cfg, tlShiftPreL_rs);title('Left, shift (baseline 100ms pre)');
cfg.channel = 'MLO11';
subplot(1,3,2);ft_singleplotER(cfg, tlShiftPreM_rs);title('Middle, shift (baseline 100ms pre)');
cfg.channel = 'MLO32';
subplot(1,3,3);ft_singleplotER(cfg, tlShiftPreR_rs);title('Right, shift (baseline 100ms pre)');
