function analysis_pilot_photo_diode(pilot)
% This function plots various aspects of the photo diode pilot. The diode
% was placed on the MEG screen at the position of the grating stimulus. The
% difference between set and actual trial length is computed, as well as
% the powerspectrum of the diode. Trigger and diode timecourses are
% plotted, together with the Hilbert transform.

% close all
if nargin<1
    pilot=15;
end
if isempty(pilot)
    pilot=15;
end


%% Load data
cfg                         = [];
if pilot==2
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode2_1200hz_20160713_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/2/pilot_diode_2.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating(1:72)+2; % in seconds (1 baseline 1 blink)
elseif pilot==3
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode3_1200hz_20160725_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/3/pilot_diode_3.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating(1:129)+2; % in seconds (1 baseline 1 blink)
elseif pilot==4
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode4_1200hz_20160725_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/4/pilot_diode_4.mat', 'log');
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==6
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode6_1200hz_20160725_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/6/pilot_diode_6.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==7
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode7_1200hz_20160726_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/7/pilot_diode_7.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==8
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode8_1200hz_20160726_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/8/pilot_diode_8.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==9
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode9_1200hz_20160726_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/9/pilot_diode_9.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==10
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode10_1200hz_20160726_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/10/pilot_diode_10.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==11
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/diodepilot11_1200hz_20160726_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/11/pilot_diode_11.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==12
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode12_1200hz_20160802_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/12/pilot_diode_12.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==13
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode13_1200hz_20160802_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/13/pilot_diode_13.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==14
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode14_1200hz_20160802_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/14/pilot_diode_14.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
elseif pilot==15
    cfg.dataset                 = '/home/electromag/matves/Data/ERF_oscillation/raw_data/pilot_photo_diode/pilotdiode15_1200hz_20160802_01.ds';
    load('/home/electromag/matves/Data/ERF_oscillation/behavioral_log/pilot_photo_diode/15/pilot_diode_15.mat', 'log');% load behavioral log file
    cfg.trialdef.poststim       = log.completeDurationGrating+2; % in seconds (1 baseline, 1 blink)
end

cfg.trialfun                = 'mytrialfun_diode'; % this is the default
cfg.trialdef.eventtype      = {'backpanel trigger', 'UPPT001'};
% the value of the stimulus trigger for photo diode pilot. These
% eventvalues were changed after this pilot.
cfg.trialdef.eventvalue     = [1, 2, 3, 4];
ifi=1/120;
cfg.trialdef.prestim        = ifi; % in seconds. Time is shifted one ifi, such that the first trigger activates at ifi instead of 0 sec.

cfg = ft_definetrial(cfg);
samples = cfg.trl;
cfg.channel = {'UPPT*', 'MEG', 'UADC001'};
% UPPT001: bitsi, UPPT002: response, UADC002: photo diode
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);
nTrials = length(data.time);

% the trialfun sees the first trigger as 0-point. For interpretation,
% display it as ifi, such that 1 sec goes from ifi:ifi:1 (instead of
% 0:ifi:1-ifi).
for iTrl = 1:nTrials
    time{iTrl} = data.time{iTrl}+ifi;
end


%{
%% Compute and plot differences setted and actual time
for i=1:length(data.time)
    diff_total(i,1) = data.time{i}(end) - (log.completeDurationGrating(1,i)+1);
    
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
cfg.channel      = 'UADC001'; % diode
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
%}

%% browse through trials
Fs=1200;

% select diode and trigger channel
for iTrl = 1:nTrials
    diode{iTrl} = data.trial{iTrl}(find(strcmp(data.label, 'UADC001')),:);
    trigger{iTrl} = data.trial{iTrl}(find(strcmp(data.label, 'UPPT001')),:);
end

% MEG Fs is 1200Hz, Screen refresh rate is 120Hz. Get the indices of the diode
% time-points in the MEG time-points (transform 1200hz to 120hz). This is
% used to shift those latencies with the executation time of the screen
% flip command
for iTrl=1:nTrials
    timeIdxDownsampled{iTrl} = (Fs*ifi)+1:Fs*ifi:length(time{iTrl});
end

% Realign data with missed Screen Flip deadlines (for illistrative purposes
% only.
triggerRealigned = trigger;
diodeRealigned = diode;

for iTrl = 1:nTrials
    idxMissedBlink = find(log.PTBtiming.blink.missedDeadline{iTrl}>0);
    idxMissedBaseline = find(log.PTBtiming.baseline.missedDeadline{iTrl}>0)+length(log.PTBtiming.blink.missedDeadline{iTrl});
    idxMissedGrating = find(log.PTBtiming.grating.missedDeadline{iTrl}>0)+length(log.PTBtiming.blink.missedDeadline{iTrl})+length(log.PTBtiming.baseline.missedDeadline{iTrl});
    idxMissed{iTrl} = [idxMissedBlink, idxMissedBaseline, idxMissedGrating]; % at which points in the data occurred a flip
    for i = 1:length(idxMissed{iTrl})
        % for every missed deadline, go to the next screen refresh latency,
        % take the rest of the data and paste it at the original flip latency
        % (where the missed deadline occurred)
        triggerRealigned{iTrl}(timeIdxDownsampled{iTrl}(idxMissed{iTrl}(i)):end-Fs*ifi) = triggerRealigned{iTrl}(timeIdxDownsampled{iTrl}(idxMissed{iTrl}(i)+1):end);
        triggerRealigned{iTrl}(timeIdxDownsampled{iTrl}(end-(i*Fs*ifi)+1:end-(i-1)*Fs*ifi))=[];
        diodeRealigned{iTrl}(timeIdxDownsampled{iTrl}(idxMissed{iTrl}(i)):end-Fs*ifi) = diodeRealigned{iTrl}(timeIdxDownsampled{iTrl}(idxMissed{iTrl}(i)+1):end);
        diodeRealigned{iTrl}(timeIdxDownsampled{iTrl}(end-(i*Fs*ifi)+1:end-(i-1)*Fs*ifi))=[];
        idxFlipAfterMissed{iTrl}(i) = idxMissed{iTrl}(i) + Fs*ifi;
    end
end
%}

%%
% add duration of Screen flip-finish to stimulus onset to the time of trigger.
timeTrigger = time;
for iTrl = 1:nTrials
    flipToStimOnsetBlink{iTrl} = log.PTBtiming.blink.StimulusOnsetTime{iTrl} - log.PTBtiming.blink.FlipTimestamp{iTrl};
    flipToStimOnsetBaseline{iTrl} = log.PTBtiming.baseline.StimulusOnsetTime{iTrl} - log.PTBtiming.baseline.FlipTimestamp{iTrl};
    flipToStimOnsetGrating{iTrl} = log.PTBtiming.grating.StimulusOnsetTime{iTrl} - log.PTBtiming.grating.FlipTimestamp{iTrl};
    flipToStimOnset{iTrl} = [flipToStimOnsetBlink{iTrl}, flipToStimOnsetBaseline{iTrl}, flipToStimOnsetGrating{iTrl}];
    % it is unclear why there are more MEG time-points than flip's. Even
    % with correction for missed flip deadlines. Crude solution: just take
    % the first n MEG timepoints.
    timeTrigger{iTrl}(timeIdxDownsampled{iTrl}(1:length(flipToStimOnset{iTrl}))) = timeTrigger{iTrl}(timeIdxDownsampled{iTrl}(1:length(flipToStimOnset{iTrl}))) + flipToStimOnset{iTrl};
end




%% Browse thorugh trials


h=figure;
% plot first trial
% plot(time{1}(1:length(diodeRealigned{1}))-ifi, diodeRealigned{1});
% NOTE: For some reason, the diode has a delay of one ifi. Possibly because
% fs_diode = ifi. So subtract one ifi from its time to synchronize with the
% trigger.
plot(time{1}-15*(1/1200),diode{1}); % align the diode and trigger data by sub-
% tracting/dividing.
hold on
% plot(timeTrigger{1}(1:length(triggerRealigned{1})), triggerRealigned{1}/10)
plot(time{1}, trigger{1}/10)
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
        plot(time{newiter}-15*(1/1200),diode{newiter});
        %         plot(time{newiter}(1:length(diodeRealigned{newiter}))-ifi,diodeRealigned{newiter});
        hold on
        plot(time{newiter}, trigger{newiter}/10)
        %         plot(timeTrigger{newiter}(1:length(triggerRealigned{newiter})), triggerRealigned{newiter}/10)
        xlim([ifi time{newiter}(end)])
        title(sprintf('trial %d', newiter))
    end
%% Calculate distribution of delays between trigger and diode

for iTrl = 1:nTrials
    % onset
    sampleOnsetTrigger(iTrl) = find(trigger{iTrl}==3, 1, 'first');
    sampleOnsetDiode(iTrl) = find(diode{iTrl}>5.7,1, 'first');
    sampleDiffOnset(iTrl) = sampleOnsetDiode(iTrl) - sampleOnsetTrigger(iTrl);
    
    % offset
    sampleOffsetTrigger(iTrl) = find(trigger{iTrl}==4, 1, 'first');
    sampleOffsetDiode(iTrl) = find(diode{iTrl}>0.8,1, 'last')+1; % take the first sample where the diode is lower than 1 again.
    sampleDiffOffset(iTrl) = sampleOffsetDiode(iTrl) - sampleOffsetTrigger(iTrl);
end
figure;
subplot(1,2,1); histogram(sampleDiffOnset); xlabel('sample diff (fs 1200hz)');title('Difference onset (diode-trigger)');
subplot(1,2,2); histogram(sampleDiffOffset); xlabel('sample diff (fs 1200hz)');title('Difference offset (diode-trigger)');

end