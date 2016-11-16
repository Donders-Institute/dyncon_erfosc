function plot_headlocation(subj, isPilot)

if nargin<1
    subj = 4;
end
if isempty(subj)
    subj = 4;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end

if isPilot
    load(sprintf('/project/3011085.02/Data/ERF_oscillation/behavioral_log/pilot/pilot_headposition_%02d.mat', subj));
else
    load(sprintf('/project/3011085.02/Data/ERF_oscillation/behavioral_log/experiment/headposition_%02d.mat', subj));
end

figure;
subplot(1,2,1)
hold on
plot(cc_dem(:,1:3)*1000) % in mm
for ii=40:40:size(cc,2)
    h = vline(ii);
    set(h, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
end
legend('x', 'y','z')
title('mean head position per trial relative to average')
xlabel('trial')
ylabel('movement (mm)')

% plot rotations
subplot(1,2,2)
hold on
plot(cc_dem(:,4:6))
for ii=40:40:size(cc,2)
    h = vline(ii);
    set(h, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
end
title('mean head orientation per trial relative to average')
xlabel('trial')
ylabel('angle (degrees)')
legend('x', 'y','z')

end