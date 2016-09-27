function [trl] = erf_osc_mytrialfun(cfg)
% This function splits up the data into trials based on trigger
% information. The triggers in the concentric grating experiment are the
% following:
%   xp.TRIG_ONSET_BLINK = 1;
%   xp.TRIG_ONSET_BASELINE = 2;
%   xp.TRIG_ONSET_GRATING = 3;
%   xp.TRIG_SHIFT = 4;
%   xp.TRIG_RESPONSE = 5;
% For analysis, grating onset is the zero time-point. The prestimulus
% baseline (trigger=2) will also be included in the trial window.
% The delay between trigger and stimulus is 15(-16) frames. To correct for
% this, select 15 frames after the trigger as the begin sample of an
% event. (event starts 15 frames after trigger).


cfg.trialdef.eventtype  = {'UPPT001'}; % define trials based on Bitsi triggers. 
cfg.trialdef.eventtyperesp = {'UPPT002'};

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'type', {'UPPT001', 'UPPT002'});
lag = 15; % in samples

% restructure events to make search easier (make a list of sample and value
% for all events)
for iEvent = 1:length(event)
    trigSample(iEvent,1) = event(iEvent).sample;
    trigValue(iEvent,1) = event(iEvent).value;
end

log = cfg.logfile.log;

trl=[];
iTrial=1;
for i = 1:length(event)
    
    isBitsiEvent = ismember(event(i).type, cfg.trialdef.eventtype);
    if isBitsiEvent
        isStartGrating = (event(i).value==3);
        if  isStartGrating
            begsample = (event(i).sample + lag) - cfg.trialdef.prestim(iTrial)*hdr.Fs;
            endsample = (event(i).sample + lag) + cfg.trialdef.poststim(iTrial)*hdr.Fs;
            offset = -cfg.trialdef.prestim(iTrial)*hdr.Fs;
            
            % get sample of baseline onset
            if ismember(event(i-1).type, cfg.trialdef.eventtype) && event(i-1).value==2 % see whether the previous event was the baseline
                sampleBaselineOnset = event(i-1).sample + lag;
            else
                idxBaselineEvent = find((trigSample<event(i).sample) & (trigValue==2), 1); % find the trigger sample with value 2 preceding grating onset
                sampleBaselineOnset = event(idxBaselineEvent).sample + lag;
            end
            
            sampleGratingOnset = event(i).sample + lag;
            % get sample of shift 
            if ismember(event(i+1).type, cfg.trialdef.eventtype) && event(i+1).value==4 % see whether the next event is the grating shift
                sampleShiftOnset = event(i+1).sample + lag;
            elseif any(iTrial==log.trlNoShift)
                sampleShiftOnset = 0;
            else
                idxShiftEvent = find((trigSample>event(i).sample) & (trigValue==4), 1); % find the next sample with value 4
                sampleShiftOnset = event(idxShiftEvent).sample + lag;
            end
            if event(i+1).sample<trigSample(end)
            isResponse = ismember(event(i+2).type, cfg.trialdef.eventtyperesp) && event(i+2).value==8;
            else
                isResponse = 0;
            end
            % find out whether there is a button response between grating onset and grating shift. 
            if find(trigSample>(sampleGratingOnset) & trigSample<(sampleShiftOnset)) % if there is a trigger in between grating onset and shift
%                 isPreviousResponse = strcmp(event(find(trigSample>(sampleGratingOnset) & trigSample<(sampleShiftOnset))).type, 'UPPT002'); % if that trigger is a button press
                    isPreviousResponse=false;
            else
                isPreviousResponse = false;
            end
            
            if event(i+1).sample<trigSample(end)
            isRespWithinTime = event(i+2).sample<endsample;
            else
                isRespWithinTime = 0;
            end
            
            
            if  isResponse && isRespWithinTime && ~isPreviousResponse
                sampleResponseOnset = event(i+2).sample + lag;
            else
                sampleResponseOnset = 0;
            end
            
           
            % some log variables
            position = log.Xpos(iTrial);

            trl(end+1,:) = [round([begsample endsample offset]) iTrial position ...
                sampleBaselineOnset sampleGratingOnset sampleShiftOnset sampleResponseOnset];
            
        
            iTrial=iTrial+1;
        end
    end
    
end
end

