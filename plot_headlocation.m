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
    load(sprintf('/project/3011085.02/clean/behavior/pilot-0%d/ses-01/pilot_headposition_0%d.mat', subj, subj));
else
    load(sprintf('/project/3011085.02/clean/behavior/subj-0%d/ses-01/headposition_0%d.mat', subj, subj));
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