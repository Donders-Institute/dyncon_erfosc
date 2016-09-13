function [trl] = mytrialfun_diode(cfg)
% first event is trial start

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'type', {'UPPT001', 'UPPT002'});

 % trialduration in samples
trl=[];
jj=1;
for i = 1:length(event)
    cfg.trialdef.eventtype  = {'UPPT001'};
    isBtsiEvent = ismember(event(i).type, cfg.trialdef.eventtype);
    if isBtsiEvent
        isStartTrial = (event(i).value==1);
        existEndTrial = jj<=length(cfg.trialdef.poststim); % is the last trial fully completed?
        if  isStartTrial && existEndTrial
            begsample = event(i).sample - cfg.trialdef.prestim*hdr.Fs;
            endsample = event(i).sample + cfg.trialdef.poststim(jj)*hdr.Fs;
            offset = -cfg.trialdef.prestim*hdr.Fs;
            grating_latency = event(i+2).sample; 
            shift_latency = event(i+3).sample;
            isResponse = strcmp(event(i+4).type, 'UPPT002');
            isRespWithinTime = event(i+4).sample<endsample;
            if  isResponse && isRespWithinTime
                resp_latency = event(i+4).sample;
            else
                resp_latency = 0;
            end
            trl(end+1,:) = [round([begsample endsample offset]) grating_latency shift_latency resp_latency];
        jj=jj+1;
        end
    end
    
end
end

