function [trl] = erf_osc_mytrialfun_sub_010(cfg)
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
for ii = 1:length(event)
    
    isBitsiEvent = ismember(event(ii).type, cfg.trialdef.eventtype);
    if isBitsiEvent
        isStartGrating = (event(ii).value==3);
        if  isStartGrating
            begsample = (event(ii).sample + lag) - cfg.trialdef.prestim(iTrial)*hdr.Fs;
            endsample = (event(ii).sample + lag) + cfg.trialdef.poststim(iTrial)*hdr.Fs;
            offset = -cfg.trialdef.prestim(iTrial)*hdr.Fs;
            
            % get sample of baseline onset
            if ismember(event(ii-1).type, cfg.trialdef.eventtype) && event(ii-1).value==2 % see whether the previous event was the baseline
                sampleBaselineOnset = event(ii-1).sample + lag;
            else
                idxBaselineEvent = find((trigSample<event(ii).sample) & (trigValue==2), 1); % find the trigger sample with value 2 preceding grating onset
                sampleBaselineOnset = event(idxBaselineEvent).sample + lag;
            end
            
            sampleGratingOnset = event(ii).sample + lag;
            % get sample of shift 
            if any(iTrial==log.trlNoShift)
                sampleShiftOnset = 0;
            elseif ismember(event(ii+1).type, cfg.trialdef.eventtype) && event(ii+1).value==4 % see whether the next event is the grating shift
                sampleShiftOnset = event(ii+1).sample + lag;
            else
                idxShiftEvent = find((trigSample>event(ii).sample) & (trigValue==4), 1); % find the next sample with value 4
                sampleShiftOnset = event(idxShiftEvent).sample + lag;
            end
            
            if ii+2<=length(event)
                isResponse = ismember(event(ii+2).type, cfg.trialdef.eventtyperesp) && event(ii+2).value==8;
                isRespWithinTime = event(ii+2).sample<endsample;
            else
                isResponse = 0;
                isRespWithinTime = 0;
            end
            
            % find out whether there is a button response between grating onset and grating shift. 
            if find(trigSample>(sampleGratingOnset) & trigSample<(sampleShiftOnset)) % if there is a trigger in between grating onset and shift
                tmpIdx = find(trigSample>(sampleGratingOnset) & trigSample<(sampleShiftOnset));
                for jj=1:length(tmpIdx)
                    isPreviousResponse = strcmp(event(tmpIdx(jj)).type, 'UPPT002'); % if that trigger is a button press
                    if isPreviousResponse
                        break
                    end
                end
            else
                isPreviousResponse = false;
            end
            
            if  isResponse && isRespWithinTime && ~isPreviousResponse
                sampleResponseOnset = event(ii+2).sample + lag;
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

