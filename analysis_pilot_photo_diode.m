function analysis_pilot_photo_diode
% This function plots various aspects of the photo diode pilot. The diode
% was placed on the MEG screen at the position of the grating stimulus. The
% difference between set and actual trial length is computed, as well as
% the powerspectrum of the diode. Trigger and diode timecourses are
% plotted, together with the Hilbert transform.


%% Load data
cfg                         = [];
cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/matspilot_1200hz_20160627_01.ds';
cfg.trialfun                = 'mytrialfun_diode'; % this is the default
cfg.trialdef.eventtype      = {'backpanel trigger', 'UPPT001'};
% the value of the stimulus trigger for photo diode pilot. These
% eventvalues were changed after this pilot.
cfg.trialdef.eventvalue     = [2, 4, 8, 16];
cfg.trialdef.prestim        = 0; % in seconds
% load behavioral log file
load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/1/test.mat', 'log');
cfg.trialdef.poststim       = log.completeDuration + 1.5; % in seconds (1.5 baseline)

cfg = ft_definetrial(cfg);
samples = cfg.trl;
cfg.channel = {'UPPT*', 'MEG', 'UADC002'};
% UPPT001: bitsi, UPPT002: response, UADC002: photo diode
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);
nTrials = length(data.time);



%% Compute and plot differences setted and actual time
for i=1:length(data.time)
    diff_total(i,1) = data.time{i}(end) - (log.completeDuration(1,i)+1.5);
    
end
figure;
subplot(1,2,1)
hist(diff_total*1000,20)
xlabel('time difference (ms)')
ylabel('#trials')
title('Time difference between measured and set time for whole trial')

%% powerspectrum diode
for ii=1:length(data.time)
    time(ii,1) = data.time{1}(end);
end

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'UADC002'; % diode
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.pad          = 6; % pad so the frequency resolution is the same for every trial (which have unequal length)
cfg.foilim          = [1 300];
powspec = ft_freqanalysis(cfg, data);
subplot(1,2,2)
plot(powspec.freq, powspec.powspctrm)
xlabel('frequency')
ylabel('power')
title('powerspectrum of the diode')

%% Display first part of all trials, timelocked to grating start

cfg=[];
cfg.trl(:,1)=data.trialinfo(:,1);
cfg.trl(:,2)=cfg.trl(:,1)+1199; % select until 1 sec after onset
cfg.trl(:,1)=data.trialinfo(:,1)-240; % select 200 ms before onset
cfg.trl(:,3)=-240; % offset
data2=ft_redefinetrial(cfg,data);

% select only trigger and diode channel
cfg=[];
cfg.channel = data2.label([1 end]);
data2=ft_selectdata(cfg,data2);

% timelock trials to onset grating.
cfg=[];
cfg.keeptrials='yes';
cfg.vartrllength=2;
tlck=ft_timelockanalysis(cfg,data2);

% plot all trials
figure;
subplot(1,3,1)
plot(tlck.time,squeeze(tlck.trial(:,2,:)))
xlabel('time (s)')
title('timelocked data diode all trials')
subplot(1,3,2)
plot(tlck.time,tlck.avg(2,:))
xlabel('time (s)')
title('timelocked data diode average all trials')

%% plot Hilbert transform
% bandpass filter around diode frequency
cfg=[];
cfg.bpfilter='yes';
cfg.bpfreq=[115 125];
cfg.bpfilttype='firws';
data3=ft_preprocessing(cfg,data2);

% hilbert transform
data4=data3;
for k = 1:numel(data4.trial)
    data4.trial{k}=abs(hilbert(data3.trial{k}'))';
end
subplot(1,3,3)
plot(data2.time{1},data3.trial{1}(2,:)) % plot bandpassed data
hold on;
plot(data2.time{1},data4.trial{1}(2,:),'r') % plot hilbert transform
xlabel('time (s)')
title('Hilbert transform diode refresh bandpass, trial 1')


%% browse through trials
nTrials = length(data.trial);
% select diode and trigger channel
for iTrl = 1:nTrials
    diode{iTrl} = data.trial{iTrl}(find(strcmp(data.label, 'UADC002')),:);
    trigger{iTrl} = data.trial{iTrl}(find(strcmp(data.label, 'UPPT001')),:);
end

h=figure;
% plot first trial
plot(data.time{1},diode{1}-0.8); % align the diode and trigger data by sub-
% tracting/dividing.
hold on
plot(data.time{1}, trigger{1}/10)
title('trial 1');
xlabel('time (s)');
newiter=[];

% create two buttons and direct a button press to the function GraphCreator
pb_next = uicontrol(h, 'Style', 'pushbutton', 'String', 'next','Position',...
    [160 0 40 40], 'CallBack', @GraphCreator, 'UserData', 1); 
pb_prev = uicontrol(h, 'Style', 'pushbutton', 'String', 'prev','Position',...
    [100 0 40 40], 'CallBack', @GraphCreator, 'UserData', 1);

    function GraphCreator(EventData, nTrials);
        % get the input of the button that has been pressed
        if pb_next.Value==1
            newiter = get(pb_next, 'UserData') +1;
            pb_next.Value=0;
        elseif pb_prev.Value==1
            newiter = get(pb_prev, 'UserData') -1;
            pb_prev.Value=0;
        end
        % replace the previous trial number by the current one in both
        % button values.
        set(pb_next, 'UserData', newiter)
        set(pb_prev, 'UserData', newiter)
        
        
        hold off
        plot(data.time{newiter},diode{newiter}-0.8);
        hold on
        plot(data.time{newiter}, trigger{newiter}/10)
        title(sprintf('trial %d', newiter))
    end

end